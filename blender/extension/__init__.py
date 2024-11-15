from pathlib import Path

import bpy
import bmesh
import mathutils

import MDAnalysis as mda

bpy.types.Scene.ProteinRunway_local_path = bpy.props.StringProperty(
    name="File",
    description="File path of the PDB with an embedded trajectory to open",
    options={"TEXTEDIT_UPDATE"},
    subtype="FILE_PATH",
    maxlen=0,
)


class ImportTrajectoryOperator(bpy.types.Operator):
    "Load the PDB"

    bl_idname = "protein_runway.import_trajectory"
    bl_label = "Import"

    # Enable undo for the operator.
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        scene = context.scene
        file_path = scene.ProteinRunway_local_path
        protein_name = Path(file_path).stem

        u = mda.Universe(file_path)

        scene = bpy.context.scene
        collection_protein = bpy.data.collections.new(protein_name)
        scene.collection.children.link(collection_protein)

        atom_mesh_object = self.draw_alpha_carbons(u, collection_protein)
        self.animate_trajectory(u, atom_mesh_object)

        # TODO (2024-11-14) This is messy and invasive, it might interfere with
        # multiple trajectories
        #
        # Fit the number of global frames to this trajectory's length.
        bpy.data.scenes[0].frame_end = len(u.trajectory)

        return {'FINISHED'}

    def draw_alpha_carbons(self, universe, collection_protein):
        atoms = universe.select_atoms('name = CA')

        center = atoms.center_of_mass()
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

        atom_mesh_object.select_set(True)

        return atom_mesh_object

    def animate_trajectory(self, universe, atom_mesh_object):
        mesh = atom_mesh_object.data

        # Snapshot first frame:
        for v in mesh.vertices:
            v.keyframe_insert('co', frame=1)
        frame = 2

        for _ in universe.trajectory[1:]:
            atoms = universe.select_atoms('name = CA')
            center = atoms.center_of_mass()
            vertices = [p - center for p in atoms.positions]

            # Perform the calculations using a BMesh:
            bm = bmesh.new()
            bm.from_mesh(mesh)
            for i, v in enumerate(bm.verts):
                v.co = mathutils.Vector(vertices[i])
            bm.to_mesh(mesh)
            bm.free()

            mesh.update()

            for v in mesh.vertices:
                v.keyframe_insert('co', frame=frame)

            frame += 1


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
        row_button.operator(ImportTrajectoryOperator.bl_idname)


def register():
    bpy.utils.register_class(ImportTrajectoryOperator)
    bpy.utils.register_class(ImportTrajectoryPanel)


def unregister():
    bpy.utils.unregister_class(ImportTrajectoryOperator)
    bpy.utils.unregister_class(ImportTrajectoryPanel)
