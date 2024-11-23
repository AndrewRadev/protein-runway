from pathlib import Path

import bpy
import bmesh

import MDAnalysis as mda

from ..lib.segmentation import generate_domain_ranges

COLORS = [
    (1.0, 1.0, 1.0, 1), # white
    (0.8, 0.1, 0.1, 1), # red
    (0.1, 0.8, 0.1, 1), # green
    (0.1, 0.1, 0.8, 1), # blue
    (0.8, 0.8, 0.1, 1), # yellow
    (0.8, 0.3, 0.8, 1), # violet
    (0.3, 0.8, 0.8, 1), # cyan
    (0.3, 0.3, 0.3, 1), # gray
]


class ImportPDBOperator(bpy.types.Operator):
    "Load the PDB"

    bl_idname = "protein_runway.import_pdb"
    bl_label  = "Import PDB"

    # Enable undo for the operator.
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        scene           = context.scene
        pdb_path        = scene.ProteinRunway_pdb_path
        add_convex_hull = scene.ProteinRunway_add_convex_hull
        protein_name    = Path(pdb_path).stem

        if len(pdb_path) == 0:
            self.report({'WARNING'}, 'No PDB provided')
            return {'CANCELLED'}

        try:
            u = mda.Universe(pdb_path)
        except ValueError as e:
            self.report({'ERROR'}, f"MDAnalysis error: {e}")
            return {'CANCELLED'}

        collection_protein = bpy.data.collections.new(protein_name)
        scene.collection.children.link(collection_protein)

        segmentations = scene.ProteinRunway_segmentations
        if len(segmentations) > 0:
            segmentation_index = scene.ProteinRunway_segmentation_index
            selected_segmentation = segmentations[segmentation_index]

            domain_regions = generate_domain_ranges(selected_segmentation.chopping)
        else:
            # One domain for the entire protein:
            domain_regions = [[range(1, len(u.atoms))]]

        atom_mesh_objects = self.draw_alpha_carbons(u, collection_protein, domain_regions, add_convex_hull)

        for o in atom_mesh_objects:
            o.select_set(True)
        bpy.ops.protein_runway.animate_trajectory('INVOKE_DEFAULT')

        # TODO (2024-11-14) This is messy and invasive, it might interfere with
        # multiple trajectories
        #
        # Fit the number of global frames to this trajectory's length.
        bpy.data.scenes[0].frame_end = len(u.trajectory)

        return {'FINISHED'}

    def draw_alpha_carbons(self, universe, collection_protein, domain_regions, add_convex_hull):
        all_atoms     = universe.select_atoms('name = CA')
        global_center = all_atoms.center_of_mass()

        radius = 0.7

        atom_mesh_objects = []

        for i, domain in enumerate(domain_regions):
            color       = COLORS[i % len(COLORS)]
            domain_name = f"Domain{i + 1:02}"

            atoms = sum(
                all_atoms.select_atoms("resnum {}:{}".format(region.start, region.stop))
                for region in domain
            )

            local_center = atoms.center_of_mass() - global_center
            vertices     = [p - global_center for p in atoms.positions]

            material = bpy.data.materials.new(f"{domain_name}_material")
            material.diffuse_color = color
            material.use_nodes = True

            # TODO: Investigate this material modification by Atomic Blender:
            #
            # mat_P_BSDF = next(n for n in material.node_tree.nodes
            #                   if n.type == "BSDF_PRINCIPLED")
            # mat_P_BSDF.inputs['Base Color'].default_value = color

            material.name = f"{domain_name}_material"

            atom_mesh = bpy.data.meshes.new(f"{domain_name}_mesh")
            atom_mesh.from_pydata(vertices, [], [])
            atom_mesh.update()
            atom_mesh_object = bpy.data.objects.new(f"{domain_name}_mesh_object", atom_mesh)

            collection_domain = bpy.data.collections.new(domain_name)
            collection_protein.children.link(collection_domain)
            collection_domain.objects.link(atom_mesh_object)

            if add_convex_hull:
                # Create the convex hull
                bm = bmesh.new()
                bm.from_mesh(atom_mesh_object.data)
                hull = bmesh.ops.convex_hull(bm, input=bm.verts)
                bm.to_mesh(atom_mesh_object.data)
                bm.free()

                # We can choose to color the hull:
                atom_mesh_object.active_material = material

                # We can make the outer points invisible:
                radius = 0

            # Create a representative UV ball for the carbons
            bpy.ops.mesh.primitive_uv_sphere_add(
                segments=32,
                ring_count=32,
                align='WORLD',
                enter_editmode=False,
                location=(0, 0, 0),
                rotation=(0, 0, 0)
            )

            representative_ball = bpy.context.view_layer.objects.active
            representative_ball.hide_set(True)
            representative_ball.scale = (radius, radius, radius)
            representative_ball.active_material = material
            representative_ball.parent = atom_mesh_object

            atom_mesh_object.instance_type = 'VERTS'

            # TODO (2024-11-14) Taken from Atomic Blender, check if it's necessary:

            # Note the collection where the ball was placed into.
            coll_all = representative_ball.users_collection
            if len(coll_all) > 0:
                coll_past = coll_all[0]
            else:
                coll_past = bpy.context.scene.collection

            # Put the atom into the new collection 'atom' and ...
            collection_protein.objects.link(representative_ball)
            # ... unlink the atom from the other collection.
            coll_past.objects.unlink(representative_ball)

            atom_mesh_objects.append(atom_mesh_object)

        return atom_mesh_objects
