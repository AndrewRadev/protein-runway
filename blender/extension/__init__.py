from pathlib import Path

import bpy
import bmesh
import mathutils

import MDAnalysis as mda

TIMER_INTERVAL = 0.001

# For PDB file input:
bpy.types.Scene.ProteinRunway_local_path = bpy.props.StringProperty(
    name="File",
    description="File path of the PDB with an embedded trajectory to open",
    options={"TEXTEDIT_UPDATE"},
    subtype="FILE_PATH",
    maxlen=0,
)

# For progress bar:
bpy.types.Scene.ProteinRunway_progress = bpy.props.FloatProperty()


class ImportPDBOperator(bpy.types.Operator):
    "Load the PDB"

    bl_idname = "protein_runway.import_pdb"
    bl_label = "Import PDB"

    # Enable undo for the operator.
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        scene        = context.scene
        file_path    = scene.ProteinRunway_local_path
        protein_name = Path(file_path).stem

        u = mda.Universe(file_path)

        collection_protein = bpy.data.collections.new(protein_name)
        scene.collection.children.link(collection_protein)

        atom_mesh_object = self.draw_alpha_carbons(u, collection_protein)

        atom_mesh_object.select_set(True)
        bpy.ops.protein_runway.animate_trajectory('INVOKE_DEFAULT')

        # TODO (2024-11-14) This is messy and invasive, it might interfere with
        # multiple trajectories
        #
        # Fit the number of global frames to this trajectory's length.
        bpy.data.scenes[0].frame_end = len(u.trajectory)

        return {'FINISHED'}

    def draw_alpha_carbons(self, universe, collection_protein):
        atoms    = universe.select_atoms('name = CA')
        center   = atoms.center_of_mass()
        vertices = [p - center for p in atoms.positions]

        color = (0.3, 0.3, 0.3, 1)
        radius = 0.7

        material = bpy.data.materials.new("Carbon_material")
        material.diffuse_color = color
        material.use_nodes = True

        # TODO: Investigate this material modification by Atomic Blender:
        #
        # mat_P_BSDF = next(n for n in material.node_tree.nodes
        #                   if n.type == "BSDF_PRINCIPLED")
        # mat_P_BSDF.inputs['Base Color'].default_value = color

        material.name = "Carbon_material"

        atom_mesh = bpy.data.meshes.new("Mesh_carbon")
        atom_mesh.from_pydata(vertices, [], [])
        atom_mesh.update()
        atom_mesh_object = bpy.data.objects.new("Carbon_mesh_object", atom_mesh)
        collection_protein.objects.link(atom_mesh_object)

        # UV balls
        bpy.ops.mesh.primitive_uv_sphere_add(
            segments=32,
            ring_count=32,
            align='WORLD',
            enter_editmode=False,
            location=(0, 0, 0),
            rotation=(0, 0, 0)
        )

        representative_ball = bpy.context.view_layer.objects.active
        representative_ball.hide_set(True)
        representative_ball.scale = (radius, radius, radius)
        representative_ball.active_material = material
        representative_ball.parent = atom_mesh_object

        atom_mesh_object.instance_type = 'VERTS'

        # TODO (2024-11-14) Taken from Atomic Blender, check if it's necessary:
        #
        # Note the collection where the ball was placed into.
        coll_all = representative_ball.users_collection
        if len(coll_all) > 0:
            coll_past = coll_all[0]
        else:
            coll_past = bpy.context.scene.collection

        # Put the atom into the new collection 'atom' and ...
        collection_protein.objects.link(representative_ball)
        # ... unlink the atom from the other collection.
        coll_past.objects.unlink(representative_ball)

        return atom_mesh_object


class AnimateTrajectoryOperator(bpy.types.Operator):
    "Import a trajectory into keyframes"

    bl_idname = "protein_runway.animate_trajectory"
    bl_label = "Animate Trajectory"

    def __init__(self):
        self.frames_inserted = 0

        self._timer = None
        self._updating = False

    def execute(self, context):
        """
        Just a necessary placeholder, the real entry point is `invoke`
        """
        return {'FINISHED'}

    def invoke(self, context, event):
        """
        Called once when the operator is invoked
        """
        atom_mesh_object = context.selected_objects[0]
        self.mesh        = atom_mesh_object.data

        pdb_path      = context.scene.ProteinRunway_local_path
        self.universe = mda.Universe(pdb_path)

        context.window_manager.modal_handler_add(self)
        self._updating = False
        self._timer = context.window_manager.event_timer_add(TIMER_INTERVAL, window=context.window)

        return {'RUNNING_MODAL'}

    def modal(self, context, event):
        """
        Main loop of the animation. Every TIMER event, a single frame is
        inserted, until they're all done. A progress bar is updated inside the
        panel.
        """
        if event.type == 'TIMER' and not self._updating:
            self._updating = True

            self.insert_frame()

            progress = self.frames_inserted / len(self.universe.trajectory)
            context.scene.ProteinRunway_progress = progress

            self._updating = False

        if self.frames_inserted >= len(self.universe.trajectory):
            context.scene.ProteinRunway_progress = progress
            return self.finish(context)

        return {'PASS_THROUGH'}

    def finish(self, context):
        context.window_manager.event_timer_remove(self._timer)
        self._timer = None
        return {'FINISHED'}

    def insert_frame(self):
        if self.frames_inserted == 0:
            # Snapshot first frame:
            for v in self.mesh.vertices:
                v.keyframe_insert('co', frame=1)
            self.frames_inserted = 1
        else:
            next(self.universe.trajectory)
            atoms = self.universe.select_atoms('name = CA')

            center   = atoms.center_of_mass()
            vertices = [p - center for p in atoms.positions]

            # Perform the calculations using a BMesh:
            bm = bmesh.new()
            bm.from_mesh(self.mesh)
            for i, v in enumerate(bm.verts):
                v.co = mathutils.Vector(vertices[i])
            bm.to_mesh(self.mesh)
            bm.free()

            self.mesh.update()

            for v in self.mesh.vertices:
                v.keyframe_insert('co', frame=(self.frames_inserted + 1))

            self.frames_inserted += 1


class ImportTrajectoryPanel(bpy.types.Panel):
    bl_label = "Import Trajectory"
    bl_idname = "PROTEINRUNWAY_PT_import_panel"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "scene"

    def draw(self, context):
        layout = self.layout
        scene = context.scene

        layout.label(text="Load a PDB File with a trajectory", icon="FILE_TICK")
        layout.separator()

        row_file = layout.row()
        row_file.prop(scene, "ProteinRunway_local_path")

        row_button = layout.row()
        row_button.operator(ImportPDBOperator.bl_idname)

        # Progress bar:
        row = self.layout.row()
        progress = context.scene.ProteinRunway_progress
        if progress == 0:
            text = "..."
        elif progress < 1:
            text = "Animating trajectory..."
        else:
            text = "Trajectory imported."
        row.progress(factor=progress, type="BAR", text=text)
        row.scale_x = 2


def register():
    bpy.utils.register_class(ImportPDBOperator)
    bpy.utils.register_class(AnimateTrajectoryOperator)
    bpy.utils.register_class(ImportTrajectoryPanel)


def unregister():
    bpy.utils.unregister_class(ImportPDBOperator)
    bpy.utils.unregister_class(AnimateTrajectoryOperator)
    bpy.utils.unregister_class(ImportTrajectoryPanel)
