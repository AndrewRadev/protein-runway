import bpy
import bmesh
import mathutils

import MDAnalysis as mda

from ..lib.segmentation import generate_domain_ranges

TIMER_INTERVAL = 0.001


class AnimateTrajectoryOperator(bpy.types.Operator):
    "Import a trajectory into keyframes"

    bl_idname = "protein_runway.animate_trajectory"
    bl_label  = "Animate Trajectory"

    def __init__(self):
        self.current_frame = 1
        self.timer         = None
        self.updating      = False

    def execute(self, context):
        """
        Just a necessary placeholder, the real entry point is `invoke`
        """
        return {'FINISHED'}

    def invoke(self, context, event):
        """
        Called once when the operator is invoked
        """
        atom_mesh_objects = context.selected_objects
        self.meshes       = [o.data for o in atom_mesh_objects]

        pdb_path      = context.scene.ProteinRunway_pdb_path
        self.universe = mda.Universe(pdb_path)

        # TODO (2024-11-19) Duplicates ImportPDBOperator
        scene = context.scene
        segmentations = scene.ProteinRunway_segmentations
        if len(segmentations) > 0:
            segmentation_index = scene.ProteinRunway_segmentation_index
            selected_segmentation = segmentations[segmentation_index]

            self.domain_regions = generate_domain_ranges(selected_segmentation.chopping)
        else:
            # One domain for the entire protein:
            self.domain_regions = [[range(1, len(u.atoms))]]

        context.window_manager.modal_handler_add(self)
        self.updating = False
        self.timer = context.window_manager.event_timer_add(TIMER_INTERVAL, window=context.window)

        return {'RUNNING_MODAL'}

    def modal(self, context, event):
        """
        Main loop of the animation. Every TIMER event, a single frame is
        inserted, until they're all done. A progress bar is updated inside the
        panel.
        """
        if event.type == 'TIMER' and not self.updating:
            self.updating = True

            self.insert_frame()

            progress = self.current_frame / len(self.universe.trajectory)
            context.scene.ProteinRunway_progress = progress

            self.updating = False

        if self.current_frame >= len(self.universe.trajectory):
            context.scene.ProteinRunway_progress = progress
            return self.finish(context)

        return {'PASS_THROUGH'}

    def finish(self, context):
        context.window_manager.event_timer_remove(self.timer)
        self.timer = None
        return {'FINISHED'}

    def insert_frame(self):
        if self.current_frame == 0:
            # Snapshot first frame:
            self.save_keyframes()
            self.current_frame += 1
        else:
            next(self.universe.trajectory)
            all_atoms     = self.universe.select_atoms('name = CA')
            global_center = all_atoms.center_of_mass()

            for i, domain in enumerate(self.domain_regions):
                atoms    = sum(all_atoms[region] for region in domain)
                mesh     = self.meshes[i]
                vertices = [p - global_center for p in atoms.positions]

                # Perform the calculations using a BMesh:
                bm = bmesh.new()
                bm.from_mesh(mesh)
                for i, v in enumerate(bm.verts):
                    v.co = mathutils.Vector(vertices[i])
                bm.to_mesh(mesh)
                bm.free()
                mesh.update()

            self.save_keyframes()
            self.current_frame += 1

    def save_keyframes(self):
        for mesh in self.meshes:
            for v in mesh.vertices:
                v.keyframe_insert('co', frame=self.current_frame)
