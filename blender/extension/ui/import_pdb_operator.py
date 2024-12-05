from pathlib import Path

import bpy
import bmesh

import MDAnalysis as mda

from .segmentations_ui_list import extract_selected_segmentation

COLORS = [
    (1.0, 1.0, 1.0, 1),  # white
    (0.8, 0.1, 0.1, 1),  # red
    (0.1, 0.8, 0.1, 1),  # green
    (0.1, 0.1, 0.8, 1),  # light blue
    (0.8, 0.8, 0.1, 1),  # yellow
    (0.8, 0.1, 0.8, 1),  # violet
    (0.1, 0.8, 0.8, 1),  # cyan
    (0.9, 0.5, 0.1, 1),  # orange
    (0.7, 0.9, 0.1, 1),  # lime
    (0.5, 0.1, 0.9, 1),  # purple
    (0.1, 0.5, 0.9, 1),  # blue
    (0.9, 0.1, 0.5, 1),  # magenta
]
GRAY = (0.3, 0.3, 0.3, 1)


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
            mda_universe = mda.Universe(pdb_path)
        except ValueError as e:
            self.report({'ERROR'}, f"MDAnalysis error: {e}")
            return {'CANCELLED'}

        self.collection_protein = bpy.data.collections.new(protein_name)
        scene.collection.children.link(self.collection_protein)

        domain_regions = extract_selected_segmentation(scene, mda_universe)

        atom_mesh_objects = self.draw_alpha_carbons(mda_universe, domain_regions, add_convex_hull)

        for o in atom_mesh_objects:
            o.select_set(True)
        bpy.ops.protein_runway.animate_trajectory('INVOKE_DEFAULT')

        # TODO (2024-11-14) This is messy and invasive, it might interfere with
        # multiple trajectories
        #
        # Fit the number of global frames to this trajectory's length.
        bpy.data.scenes[0].frame_end = len(mda_universe.trajectory)

        return {'FINISHED'}

    def draw_alpha_carbons(self, universe, domain_regions, add_convex_hull):
        all_atoms     = universe.select_atoms('name = CA')
        global_center = all_atoms.center_of_mass()
        all_resnums   = {a.resnum for a in all_atoms.atoms}

        atom_mesh_objects = []
        seen_resnums = set()

        for i, domain in enumerate(domain_regions):
            color       = COLORS[i % len(COLORS)]
            domain_name = f"Domain{i + 1:02}"

            atom_group = sum(
                all_atoms.select_atoms("resnum {}:{}".format(region.start, region.stop))
                for region in domain
            )
            seen_resnums |= {a.resnum for a in atom_group.atoms}
            vertices = [p - global_center for p in atom_group.positions]

            atom_mesh_object = self.draw_domain(vertices, domain_name, color, add_convex_hull)
            atom_mesh_objects.append(atom_mesh_object)

        # Add leftover atoms, if any:
        unseen_resnums = all_resnums - seen_resnums
        if len(unseen_resnums) > 0:
            color       = GRAY
            domain_name = "Domain_NA"

            atom_group = sum(
                all_atoms.select_atoms("resnum {}".format(str(resnum)))
                for resnum in unseen_resnums
            )
            vertices = [p - global_center for p in atom_group.positions]

            # Note: We do not add a convex hull for the unclassified atoms,
            # since they will be spread out
            atom_mesh_object = self.draw_domain(vertices, domain_name, color, False)
            atom_mesh_objects.append(atom_mesh_object)

        return atom_mesh_objects

    def draw_domain(self, vertices, domain_name, color, add_convex_hull):
        material = bpy.data.materials.new(f"{domain_name}_material")
        material.diffuse_color = color
        material.use_nodes = True

        material.name = f"{domain_name}_material"

        atom_mesh = bpy.data.meshes.new(f"{domain_name}_mesh")
        atom_mesh.from_pydata(vertices, [], [])
        atom_mesh.update()
        atom_mesh_object = bpy.data.objects.new(f"{domain_name}_mesh_object", atom_mesh)

        collection_domain = bpy.data.collections.new(domain_name)
        self.collection_protein.children.link(collection_domain)
        collection_domain.objects.link(atom_mesh_object)

        if add_convex_hull:
            # Create the convex hull
            bm = bmesh.new()
            bm.from_mesh(atom_mesh_object.data)
            bmesh.ops.convex_hull(bm, input=bm.verts)
            bm.to_mesh(atom_mesh_object.data)
            bm.free()

            # We can choose to color the hull:
            atom_mesh_object.active_material = material

            # We can make the outer points invisible:
            radius = 0.0
        else:
            radius = 1.0

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

        self.collection_protein.objects.link(representative_ball)

        return atom_mesh_object
