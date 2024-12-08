import bpy
import bmesh
import mathutils

import MDAnalysis as mda

from .segmentations_ui_list import extract_selected_segmentation

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
        scene             = context.scene
        atom_mesh_objects = context.selected_objects
        self.meshes       = [o.data for o in atom_mesh_objects]

        pdb_path            = context.scene.ProteinRunway_pdb_path
        self.universe       = mda.Universe(pdb_path)
        self.domain_regions = extract_selected_segmentation(scene, self.universe)

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
            all_resnums   = {a.resnum for a in all_atoms.atoms}
            global_center = all_atoms.center_of_mass()
            seen_resnums  = set()

            for i, domain in enumerate(self.domain_regions):
                atom_group = sum(
                    all_atoms.select_atoms("resnum {}:{}".format(region.start, region.stop))
                    for region in domain
                )
                seen_resnums |= {a.resnum for a in atom_group.atoms}

                mesh     = self.meshes[i]
                vertices = [p - global_center for p in atom_group.positions]

                # Perform the calculations using a BMesh:
                self.update_mesh_coordinates(mesh, vertices)

            # Add leftover atoms, if any:
            unseen_resnums = all_resnums - seen_resnums
            if len(unseen_resnums) > 0:
                atom_group = sum(
                    all_atoms.select_atoms("resnum {}".format(str(resnum)))
                    for resnum in unseen_resnums
                )

                mesh     = self.meshes[len(self.domain_regions)]
                vertices = [p - global_center for p in atom_group.positions]

                self.update_mesh_coordinates(mesh, vertices)

            self.save_keyframes()
            self.current_frame += 1

    def update_mesh_coordinates(self, mesh, vertices):
        # Perform the calculations using a BMesh:
        bm = bmesh.new()
        bm.from_mesh(mesh)
        for i, v in enumerate(bm.verts):
            v.co = mathutils.Vector(vertices[i])

        # Write the bmesh to the real mesh
        bm.to_mesh(mesh)
        bm.free()
        mesh.update()

    def save_keyframes(self):
        for mesh in self.meshes:
            for v in mesh.vertices:
                v.keyframe_insert('co', frame=self.current_frame)
