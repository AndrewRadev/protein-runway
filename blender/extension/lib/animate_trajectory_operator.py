import bpy
import bmesh
import mathutils

import MDAnalysis as mda

TIMER_INTERVAL = 0.001


class AnimateTrajectoryOperator(bpy.types.Operator):
    "Import a trajectory into keyframes"

    bl_idname = "protein_runway.animate_trajectory"
    bl_label  = "Animate Trajectory"

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
