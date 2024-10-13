import bpy
import bmesh
import math
import logging


class CreateWaterOperator(bpy.types.Operator):
    bl_idname = "protein_runway.create_water_operator"
    bl_label = "Create a water molecule"

    def execute(self, context):
        red  = (1, 0, 0, 1)
        gray = (0.3, 0.3, 0.3, 1)

        o  = self.create_atom("Oxygen",     red,  1.0)
        h1 = self.create_atom("Hydrogen_1", gray, 0.7)
        h2 = self.create_atom("Hydrogen_2", gray, 0.7)
        # b1 = self.create_bond(o, h1, 0.1)
        # b2 = self.create_bond(o, h2, 0.1)

        self.set_atom_position(o,  (0.0, 0.0,  0.5),  frame=1)
        self.set_atom_position(h1, (0.0, 1.5,  -1.0), frame=1)
        self.set_atom_position(h2, (0.0, -1.5, -1.0), frame=1)
        # self.update_bond(b1, o, h1, frame=1)
        # self.update_bond(b2, o, h2, frame=1)

        self.set_atom_position(o,  (0.0, 0.0,  1.5),  frame=10)
        self.set_atom_position(h1, (0.0, 1.5,  -0.0), frame=10)
        self.set_atom_position(h2, (0.0, -1.5, -0.0), frame=10)
        # self.update_bond(b1, o, h1, frame=10)
        # self.update_bond(b2, o, h2, frame=10)

        self.set_atom_position(o,  (0.0, 0.0,  0.5),  frame=20)
        self.set_atom_position(h1, (0.0, 1.5,  -1.0), frame=20)
        self.set_atom_position(h2, (0.0, -1.5, -1.0), frame=20)
        # self.update_bond(b1, o, h1, frame=20)
        # self.update_bond(b2, o, h2, frame=20)

        for item in [o, h1, h2]:
            item.select_set(True)

        return {'FINISHED'}

    def create_atom(self, name, color, radius):
        # Create an empty mesh and the object.
        mesh = bpy.data.meshes.new(name)
        atom = bpy.data.objects.new(name, mesh)

        material = bpy.data.materials.new(f"{name}_material")

        material.diffuse_color = color
        atom.active_material   = material

        # Add the object into the scene.
        bpy.context.collection.objects.link(atom)

        # Select the newly created object
        bpy.context.view_layer.objects.active = atom
        atom.select_set(True)

        # Construct the bmesh sphere and assign it to the blender mesh.
        bm = bmesh.new()
        bmesh.ops.create_uvsphere(bm, u_segments=32, v_segments=16, radius=radius)
        bm.to_mesh(mesh)
        bm.free()

        bpy.ops.object.modifier_add(type='SUBSURF')
        bpy.ops.object.shade_smooth()

        return atom

    def create_bond(self, object1, object2, radius):
        vector = (object2.location - object1.location)
        distance = vector.length

        cylinder = bpy.ops.mesh.primitive_cylinder_add(
            radius=radius,
            depth=distance,
        )
        bpy.context.active_object.name = f"Bond {object1.name}-{object2.name}"

        return cylinder

    def set_atom_position(self, atom, location, frame):
        atom.location = location
        atom.keyframe_insert(data_path="location", frame=frame)

    def update_bond(self, bond, object1, object2, frame):
        vector = (object2.location - object1.location)
        distance = vector.length

        phi = math.atan2(vector[1], vector[0])
        theta = math.acos(vector[2]/distance)

        # Depth of the cylinder:
        # bond.dimensions[2] = distance

        bond.location = object1.location + vector / 2
        bond.rotation = (0, theta, phi)

        # bond.keyframe_insert(data_path="dimensions", frame=frame)
        bond.keyframe_insert(data_path="location", frame=frame)
        bond.keyframe_insert(data_path="rotation", frame=frame)




def draw_create_water_menu(self, context):
    self.layout.operator(CreateWaterOperator.bl_idname, text="Create Water 💧")


def register():
    bpy.utils.register_class(CreateWaterOperator)


def unregister():
    bpy.utils.unregister_class(CreateWaterOperator)


if __name__ == "__main__":
    register()

    bpy.types.TOPBAR_MT_edit.append(draw_create_water_menu)

    bpy.data.scenes[0].frame_end = 20
