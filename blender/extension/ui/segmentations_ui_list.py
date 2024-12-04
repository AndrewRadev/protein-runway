import bpy

from ..lib.segmentation import generate_domain_ranges

class SegmentationMethodItem(bpy.types.PropertyGroup):
    """
    A single segmentation method, wrapped in a custom class
    """
    name: bpy.props.StringProperty(name="Method name")


class SegmentationItem(bpy.types.PropertyGroup):
    """
    Group of properties representing an item in the list of domain counts to
    pick from for a given method.
    """
    method_name:  bpy.props.StringProperty(name="Method name")
    domain_count: bpy.props.StringProperty(name="Domain count")
    chopping:     bpy.props.StringProperty(name="Chopping")


class SegmentationMethodsUiList(bpy.types.UIList):
    """
    List of unique segmentation methods
    """

    bl_idname = "PROTEINRUNWAY_UL_segmentation_methods_ui_list"

    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        layout.label(text=item.name, icon='GRAPH')


class SegmentationParamsUiList(bpy.types.UIList):
    """
    List of possible segmentations of domains, represented by the domain count,
    but bound to a particular method.
    """

    bl_idname = "PROTEINRUNWAY_UL_segmentation_params_ui_list"

    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        layout.label(text=item.domain_count, icon='OPTIONS')

    def filter_items(self, context, data, propname):
        scene = context.scene

        active_method_index = scene.ProteinRunway_segmentation_method_index
        active_method       = scene.ProteinRunway_segmentation_methods[active_method_index]

        active_params_index = scene.ProteinRunway_segmentation_params_index
        active_item         = scene.ProteinRunway_segmentation_items[active_params_index]

        items = scene.ProteinRunway_segmentation_items
        filter_flags = [self.bitflag_filter_item] * len(items)

        for index, item in enumerate(items):
            if active_method.name != item.method_name:
                filter_flags[index] &= True

        return filter_flags, []


def extract_selected_segmentation(scene, mda_universe):
    active_method_index       = scene.ProteinRunway_segmentation_method_index
    active_segmentation_index = scene.ProteinRunway_segmentation_params_index

    if len(scene.ProteinRunway_segmentation_items) > 0:
        active_segmentation = scene.ProteinRunway_segmentation_items[active_segmentation_index]
        active_method       = scene.ProteinRunway_segmentation_methods[active_method_index]

        if active_method.name != active_segmentation.method_name:
            # then the "selected" one in the list is actually hidden, so let's
            # just take the first one that is relevant to this method:
            active_segmentation = next((
                item
                for item in scene.ProteinRunway_segmentation_items
                if item.method_name == active_method.name
            ), None)
    else:
        active_segmentation = None

    if active_segmentation is not None and active_segmentation.chopping != '':
        domain_regions = generate_domain_ranges(active_segmentation.chopping)
    else:
        # One domain for the entire protein:
        domain_regions = [[range(1, max(a.resnum for a in mda_universe.atoms))]]

    return domain_regions
