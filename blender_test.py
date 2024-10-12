import bpy

class HelloWorldPanel(bpy.types.Panel):
    """
    Creates a Panel in the Object properties window
    """
    bl_label       = "Hello World Panel"
    bl_idname      = "protein_runway_hello_panel"
    bl_space_type  = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context     = "object"

    def draw(self, context):
        layout = self.layout
        obj = context.object

        # Hello world text
        row = layout.row()
        row.label(text="Hello world!", icon='WORLD_DATA')

        # Object name read/write
        row = layout.row()
        row.label(text="Active object is: " + obj.name)
        row = layout.row()
        row.prop(obj, "name")

        # Add cube
        row = layout.row()
        row.operator("mesh.primitive_cube_add")

        # Move up on the Z axis:
        row = layout.row()
        row.operator(MoveUpOperator.bl_idname)


class MoveUpOperator(bpy.types.Operator):
    bl_idname = f"protein_runway.move_up_operator"
    bl_label = "Move Up (on the Z-axis)"

    def execute(self, context):
        bpy.context.object.location.z += 2.0
        self.report({'INFO'}, f"Moving the object up on the z-axis")
        return {'FINISHED'}


def register():
    bpy.utils.register_class(HelloWorldPanel)
    bpy.utils.register_class(MoveUpOperator)


def unregister():
    bpy.utils.unregister_class(HelloWorldPanel)
    bpy.utils.unregister_class(MoveUpOperator)


if __name__ == "__main__":
    register()
