import bpy

class SegmentationItem(bpy.types.PropertyGroup):
    """
    Group of properties representing an item in the segmentation list
    """
    method: bpy.props.StringProperty(
        name="Method",
        description="Segmentation method",
        default="",
    )

    chopping: bpy.props.StringProperty(
        name="Chopping",
        description="Residue ranges grouped in domains",
        default="",
    )

class SegmentationsUiList(bpy.types.UIList):
    """
    List of possible segmentations of domains
    """

    bl_idname = "PROTEINRUNWAY_UL_segmentations_ui_list"

    def draw_item(
        self,
        context,
        layout,
        data,
        item,
        icon,
        active_data,
        active_propname,
        index,
    ):
        layout.label(text=item.method, icon='FILE_3D')
