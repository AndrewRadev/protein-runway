import bpy

class SegmentationMethodItem(bpy.types.PropertyGroup):
    """
    Group of properties representing an item in the segmentatio methods list
    """
    method: bpy.props.StringProperty(
        name="Method",
        description="Segmentation method",
    )

    domain_counts: bpy.props.StringProperty(
        name="DomainCounts",
        description="Number of domains for a given method, comma-separated",
    )

class SegmentationDomainCountItem(bpy.types.PropertyGroup):
    """
    Group of properties representing an item in the list of domain counts to
    pick from for a given method.
    """
    domain_count: bpy.props.StringProperty(
        name="DomainCount",
        description="Number of domains for a given method",
    )

    chopping: bpy.props.StringProperty(
        name="Chopping",
        description="Residue ranges grouped in domains",
    )

class SegmentationMethodsUiList(bpy.types.UIList):
    """
    List of possible segmentations of domains
    """

    bl_idname = "PROTEINRUNWAY_UL_segmentation_methods_ui_list"

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

class SegmentationDomainCountsUiList(bpy.types.UIList):
    """
    List of possible segmentations of domains
    """

    bl_idname = "PROTEINRUNWAY_UL_segmentation_domain_counts_ui_list"

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
        layout.label(text=item.domain_count, icon='FILE_3D')
