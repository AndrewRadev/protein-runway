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
        b1 = self.create_bond(o, h1, 0.1)
        b2 = self.create_bond(o, h2, 0.1)

        # Starting position
        self.set_atom_position(o,  (0.0, 0.0,  0.5),  frame=1)
        self.set_atom_position(h1, (0.0, 1.5,  -1.0), frame=1)
        self.set_atom_position(h2, (0.0, -1.5, -1.0), frame=1)

        # Up the z-axis
        self.set_atom_position(o,  (0.0, 0.0,  1.5),  frame=10)
        self.set_atom_position(h1, (0.0, 1.5,  -0.0), frame=10)
        self.set_atom_position(h2, (0.0, -1.5, -0.0), frame=10)

        # Stretch the hydrogens apart
        self.set_atom_position(o,  (0.0, 0.0,  0.5),  frame=20)
        self.set_atom_position(h1, (0.0, 2.0,  -0.7), frame=20)
        self.set_atom_position(h2, (0.0, -2.0, -0.7), frame=20)

        # Put the hydrogens back together
        self.set_atom_position(o,  (0.0, 0.0,  1.5),  frame=30)
        self.set_atom_position(h1, (0.0, 1.5,  -0.0), frame=30)
        self.set_atom_position(h2, (0.0, -1.5, -0.0), frame=30)

        # Down the z-axis
        self.set_atom_position(o,  (0.0, 0.0,  0.5),  frame=40)
        self.set_atom_position(h1, (0.0, 1.5,  -1.0), frame=40)
        self.set_atom_position(h2, (0.0, -1.5, -1.0), frame=40)

        for item in [o, h1, h2, b1, b2]:
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

    def create_bond(self, atom1, atom2, radius):
        # Create new connector mesh and mesh object and link to scene
        m = bpy.data.meshes.new('connector')

        bm = bmesh.new()
        v1 = bm.verts.new(atom1.location)
        v2 = bm.verts.new(atom2.location)
        e  = bm.edges.new([v1, v2])

        bm.to_mesh(m)

        name = f"Bond {atom1.name}-{atom2.name}"
        bond = bpy.data.objects.new(name, m)
        bpy.context.scene.collection.objects.link(bond)

        # Hook connector vertices to respective atoms
        for i, atom in enumerate([atom1, atom2]):
            bpy.ops.object.select_all(action = 'DESELECT')
            atom.select_set(True)
            bond.select_set(True)
            bpy.context.view_layer.objects.active = bond # Set connector as active

            # Select vertex
            bpy.ops.object.mode_set(mode='OBJECT')
            bond.data.vertices[i].select = True
            bpy.ops.object.mode_set(mode='EDIT')

            bpy.ops.object.hook_add_selob() # Hook to cylinder

            bpy.ops.object.mode_set(mode='OBJECT')
            bond.data.vertices[i].select = False

        m = bond.modifiers.new('Skin', 'SKIN')

        # Set radius:
        for v in bond.data.skin_vertices[0].data:
            v.radius = (radius, radius)

        # Add shading:
        m.use_smooth_shade = True
        m = bond.modifiers.new('Subsurf', 'SUBSURF')
        m.levels = 2
        m.render_levels = 2

        return bond

    def set_atom_position(self, atom, location, frame):
        atom.location = location
        atom.keyframe_insert(data_path="location", frame=frame)


def draw_create_water_menu(self, context):
    self.layout.operator(CreateWaterOperator.bl_idname, text="Create Water 💧")


def register():
    bpy.utils.register_class(CreateWaterOperator)


def unregister():
    bpy.utils.unregister_class(CreateWaterOperator)


if __name__ == "__main__":
    register()

    bpy.types.TOPBAR_MT_edit.append(draw_create_water_menu)

    bpy.data.scenes[0].frame_end = 40
