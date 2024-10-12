import bpy
import bmesh
import math
import logging


class CreateWaterOperator(bpy.types.Operator):
    bl_idname = "protein_runway.create_water_operator"
    bl_label = "Create a water molecule"

    def create_atom(self, name, color, radius, location):
        # Create an empty mesh and the object.
        mesh = bpy.data.meshes.new(name)
        atom = bpy.data.objects.new(name, mesh)

        material = bpy.data.materials.new(f"{name}_material")

        material.diffuse_color = color
        atom.active_material   = material
        atom.location          = location

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

        phi = math.atan2(vector[1], vector[0])
        theta = math.acos(vector[2]/distance)

        cylinder = bpy.ops.mesh.primitive_cylinder_add(
            radius = radius,
            depth = distance,
            location = object1.location + vector / 2,
            rotation = (0, theta, phi)
        )

        bpy.context.active_object.name = f"Bond {object1.name}-{object2.name}"

        return cylinder


    def execute(self, context):
        red  = (1, 0, 0, 1)
        gray = (0.3, 0.3, 0.3, 1)

        o  = self.create_atom("Oxygen",     red,  1.0, (0.0, 0.0,  0.5))
        h1 = self.create_atom("Hydrogen_1", gray, 0.7, (0.0, 1.5,  -1.0))
        h2 = self.create_atom("Hydrogen_2", gray, 0.7, (0.0, -1.5, -1.0))

        b1 = self.create_bond(o, h1, 0.1)
        b2 = self.create_bond(o, h2, 0.1)

        for item in [o, h1, h2, b1, b2]:
            item.select_set(True)

        return {'FINISHED'}


def draw_create_water_menu(self, context):
    self.layout.operator(CreateWaterOperator.bl_idname, text="Create Water 💧")


def register():
    bpy.utils.register_class(CreateWaterOperator)


def unregister():
    bpy.utils.unregister_class(CreateWaterOperator)


if __name__ == "__main__":
    register()
    bpy.types.TOPBAR_MT_edit.append(draw_create_water_menu)
