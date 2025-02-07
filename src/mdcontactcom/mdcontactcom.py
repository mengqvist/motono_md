import argparse
import os
import sys

from argparse import RawTextHelpFormatter
from contact_profile_calculator import profile_main, ALL_ATOM_CONTACT_CUTOFF,\
                                               ONLY_HEAVY_ATOM_CONTACT_CUTOFF,\
                                               CA_ATOM_CONTACT_CUTOFF, CB_ATOM_CONTACT_CUTOFF
from contact_similarity_calculator import similarity_main
from similarity_drawer import drawer_main, SIMILARITY_TANIMOTO, SIMILARITY_EUCLIDEAN
from similarity_plotter import main as plotter_main
from os.path import join, basename, isfile, exists


def run(args):
    """system execution function

    Args:
        args (NameSpace): Argument parsed namespace
    """
    trajectory1_pdb = args.trajectory1
    trajectory2_pdb = args.trajectory2
    trajectory1_name = args.traj1_name
    trajectory2_name = args.traj2_name
    region = args.region
    contact_cutoff = None
    num_process = args.num_process
    out_dist_map = args.out_dist_map
    output_dir = args.output_dir
    contact_file1 = join(output_dir, f'{trajectory1_name}_contact_profile.csv')
    contact_file2 = join(output_dir, f'{trajectory2_name}_contact_profile.csv')
    distance_file1 = join(output_dir, f'{trajectory1_name}_distance_profile.csv')
    distance_file2 = join(output_dir, f'{trajectory2_name}_distance_profile.csv')
    log_file_path1 = os.path.join(output_dir, f'{trajectory1_name}_contact_profile_calculator.log')
    log_file_path2 = os.path.join(output_dir, f'{trajectory2_name}_contact_profile_calculator.log')

    # calculate the contact frequencies
    for trajectory, contact_file, dist_file, log_file in zip([trajectory1_pdb, trajectory2_pdb],
                                   [contact_file1, contact_file2],
                                   [distance_file1, distance_file2],
                                   [log_file_path1, log_file_path2]):

        profile_main(trajectory=trajectory,
                     region=region,
                     contact_cutoff=contact_cutoff,
                     num_process=num_process,
                     out_dist_map=out_dist_map,
                     target_atom=args.target_atom,
                     contact_file=contact_file,
                     distance_file=dist_file,
                     log_file_path=log_file)

    # calculate similarity
    similarity_file = join(output_dir, f'{trajectory1_name}_vs_{trajectory2_name}_similarity_table.csv')
    diff_file = join(output_dir, f'{trajectory1_name}_vs_{trajectory2_name}_contact_profile_diff.csv')
    log_file = join(output_dir, f'{trajectory1_name}_vs_{trajectory2_name}_contact_similarity_calculator.log')

    similarity_main(cp1=contact_file1,
                    cp2=contact_file2,
                    diff_file=diff_file,
                    similarity_file=similarity_file,
                    log_file=log_file)

    # drawer_main(
    #     contact_profile_diff_file=diff_file,
    #     similarity_table_file=similarity_file,
    #     trajectory=trajectory1_pdb,
    #     exclude_region='',
    #     similarity_mode=args.similarity,
    #     output_dir=os.path.join(output_dir, args.output_dir))

    # plotter_main(
    #     similarity_table_file=os.path.join(output_dir, 'similarity_table.csv'),
    #     bin=50,
    #     exclude_region='',
    #     output_dir = output_dir)

    #
    # plotter_main(
    #     similarity_table_file=args.similarity_table,
    #     bin=args.bin,
    #     exclude_region=args.noregion,
    #     output_dir=args.output_dir)


# def profile(args):
#     """contact_profile_calculator execution function
#
#     Args:
#         args (NameSpace): Argument parsed namespace
#     """
#     trajectory = args.trajectory
#     topology = args.topology
#     contact_cutoff = None
#     output_dir = os.path.splitext(os.path.basename(trajectory))[0]
#
#     # Trajectory conversion processing
#     trajectory_pdb = trajectory
#
#     profile_main(trajectory=trajectory_pdb,
#                  region=args.region,
#                  contact_cutoff=contact_cutoff,
#                  num_process=args.num_process,
#                  out_dist_map=args.out_dist_map,
#                  target_atom=args.target_atom,
#                  output_dir=output_dir)



def get_argparser():
    parser = argparse.ArgumentParser(
        description='MDContactCom System'
    )

    subparsers = parser.add_subparsers(
        dest='subparser_name'
    )

    # system execution command
    parser_run = subparsers.add_parser(
        'run',
        description="Run MDContactCom Systems",
        formatter_class=RawTextHelpFormatter,
        help='see `run -h`'
    )
    parser_run.add_argument(
        'trajectory1',
        help='trajectory file is loaded by filename extension.\n'
             'Supported trajectory formats include ".pdb".'
    )
    parser_run.add_argument(
        'trajectory2',
        help='trajectory file is loaded by filename extension.\n'
             'Supported trajectory formats include ".pdb".'
    )
    parser_run.add_argument(
        'traj1_name',
        help='a short string describing the trajectory 1 file, used for downstream file naming'
    )
    parser_run.add_argument(
        'traj2_name',
        help='a short string describing the trajectory 1 file, used for downstream file naming'
    )
    parser_run.add_argument(
        '-r', '--region',
        help='Contact profile calculation target area\n'
             'Format: (<chain>):(<resnum>) \n'
             'Multiple specifications can be specified by separating them with ","',
        default=''
    )
    parser_run.add_argument(
        '-c', '--contact_cutoff',
        help='Contact judgment threshold.\n'
             'Option "-t 0". Defaults to  {0}.\n'
             'Option "-t 1". Defaults to  {1}.\n'
             'Option "-t 2". Defaults to  {2}.\n'
             'Option "-t 3". Defaults to  {3}.\n'.format(
                 ALL_ATOM_CONTACT_CUTOFF,
                 ONLY_HEAVY_ATOM_CONTACT_CUTOFF,
                 CA_ATOM_CONTACT_CUTOFF,
                 CB_ATOM_CONTACT_CUTOFF
             ),
        type=float
    )
    parser_run.add_argument(
        '-p', '--num_process',
        help='Number of processes.(Run on multiprocessing)',
        default=1,
        type=int
    )
    parser_run.add_argument(
        '-m', '--out_dist_map',
        help='Output an intermediate file (distance matrix).',
        action='store_true'
    )
    parser_run.add_argument(
        '-t', '--target_atom',
        help='Interaction analysis target atoms.\n'
             '0: All atoms\n'
             '1: Heavy atoms (Atoms excluding hydrogen atoms)[Default]\n'
             '2: CA atoms\n'
             '3: CB atoms\n',
        default=1,
        choices=[0, 1, 2, 3],
        type=int
    )
    parser_run.add_argument(
        '-s', '--similarity',
        help='Similarity coefficient to be visualized',
        type=str,
        choices=[SIMILARITY_TANIMOTO, SIMILARITY_EUCLIDEAN],
        default=SIMILARITY_TANIMOTO
    )
    parser_run.add_argument(
        '-o', '--output_dir',
        help='Directory in which to save output files',
        default='outputs'
    )
    parser_run.set_defaults(handler=run)

    # contact_profile_calculator execution command
    parser_profile = subparsers.add_parser(
        'profile',
        description="Run contact_profile_calculator",
        formatter_class=RawTextHelpFormatter,
        help='see `profile -h`'
    )
    parser_profile.add_argument(
        'trajectory',
        help='trajectory file is loaded by filename extension.\n'
             'Supported trajectory formats include ".pdb", ".dcd", ".dtr", ".xtc", ".mdcrd".'
    )
    parser_profile.add_argument(
        '-r', '--region',
        help='Contact profile calculation target area\n'
             'Format: (<chain>):(<resnum>) \n'
             'Multiple specifications can be specified by separating them with ","',
        default=''
    )
    parser_profile.add_argument(
        '-c', '--contact_cutoff',
        help='Contact judgment threshold.\n'
             'Option "-t 0". Defaults to  {0}.\n'
             'Option "-t 1". Defaults to  {1}.\n'
             'Option "-t 2". Defaults to  {2}.\n'
             'Option "-t 3". Defaults to  {3}.\n'.format(
                 ALL_ATOM_CONTACT_CUTOFF,
                 ONLY_HEAVY_ATOM_CONTACT_CUTOFF,
                 CA_ATOM_CONTACT_CUTOFF,
                 CB_ATOM_CONTACT_CUTOFF
             ),
        type=float
    )
    parser_profile.add_argument(
        '-p', '--num_process',
        help='Number of processes.(Run on multiprocessing)',
        default=1,
        type=int
    )
    parser_profile.add_argument(
        '-m', '--out_dist_map',
        help='Output an intermediate file (distance matrix).',
        action='store_true'
    )
    parser_profile.add_argument(
        '-t', '--target_atom',
        help='Interaction analysis target atoms.\n'
             '0: All atoms\n'
             '1: Heavy atoms (Atoms excluding hydrogen atoms)[Default]\n'
             '2: CA atoms\n'
             '3: CB atoms\n',
        default=1,
        choices=[0, 1, 2, 3],
        type=int
    )
    parser_profile.set_defaults(handler=run)

    # # contact_similarity_calculator
    # parser_similarity = subparsers.add_parser(
    #     'similarity',
    #     description="Run contact_similarity_calculator",
    #     formatter_class=RawTextHelpFormatter,
    #     help='see `similarity -h`'
    # )
    # parser_similarity.add_argument(
    #     'contact_profile_1',
    #     help='contact profile file',
    #     type=str
    # )
    # parser_similarity.add_argument(
    #     'contact_profile_2',
    #     help='contact profile file',
    #     type=str
    # )
    # parser_similarity.set_defaults(handler=similarity)

    # # similarity_drawer
    # parser_drawer = subparsers.add_parser(
    #     'drawer',
    #     description="Run similarity_drawer",
    #     formatter_class=RawTextHelpFormatter,
    #     help='see `drawer -h`'
    # )
    # parser_drawer.add_argument(
    #     'contact_profile',
    #     help='contact profile diff file.'
    # )
    # parser_drawer.add_argument(
    #     'similarity_table',
    #     help='Similarity table file.'
    # )
    # parser_drawer.add_argument(
    #     'trajectory',
    #     help='trajectory file is loaded by filename extension.\n'
    #          'Supported trajectory formats include ".pdb", ".dcd", ".dtr", ".xtc", ".mdcrd".'
    # )
    # parser_drawer.add_argument(
    #     '--topology', '-top',
    #     help='topology file is loaded by filename extension.\n'
    #          'Supported topology formats include ".pdb", ".cms", ".psf", ".prmtop".',
    #     type=str,
    #     default=None
    # )
    # parser_drawer.add_argument(
    #     '--noregion', '-n',
    #     help='Exclude residues.\n'
    #          'Format: (<chain>):(<resnum>) \n'
    #          'Multiple specifications can be specified by separating them with ","',
    #     default=''
    # )
    # parser_drawer.add_argument(
    #     '-s', '--similarity',
    #     help='Similarity coefficient to be visualized',
    #     type=str,
    #     choices=[SIMILARITY_TANIMOTO, SIMILARITY_EUCLIDEAN],
    #     default=SIMILARITY_TANIMOTO
    # )
    # parser_drawer.set_defaults(handler=drawer)

    # # similarity_plotter
    # parser_plotter = subparsers.add_parser(
    #     'plotter',
    #     description="Run similarity_plotter",
    #     formatter_class=RawTextHelpFormatter,
    #     help='see `plotter -h`'
    # )
    # parser_plotter.add_argument(
    #     'similarity_table',
    #     help='similarity table file',
    #     type=str
    # )
    # parser_plotter.add_argument(
    #     '--bin', '-b',
    #     help='Display width of graph scale. (Defalt 50)',
    #     type=int,
    #     default=50
    # )
    # parser_plotter.add_argument(
    #     '--noregion', '-n',
    #     help='Exclude residues.\n'
    #          'Format: (<chain>):(<resnum>) \n'
    #          'Multiple specifications can be specified by separating them with ","',
    #     default=''
    # )
    # parser_plotter.set_defaults(handler=plotter)

    return parser


def not_in_sub_cmd(argv: list) -> bool:
    """A function that determines whether a subcommand is specified in the standard input

    Args:
        argv (list): Standard input when executing a command: sys.argv

    Returns:
        bool: If the specified subcommand is not included TRUE
    """
    if(len(argv) < 2):
        return True

    subcommands = ['-h', 'run', 'profile', 'similarity', 'drawer', 'plotter']
    return argv[1] not in subcommands


if __name__ == '__main__':
    parser = get_argparser()
    argv = None
    if(not_in_sub_cmd(sys.argv)):
        argv = list(['run'])
        argv.extend(sys.argv[1:])
    else:
        argv = sys.argv[1:]

    args_parser = parser.parse_args(argv)
    if hasattr(args_parser, 'handler'):
        # Execution of the set handler function
        args_parser.handler(args_parser)
    else:
        # Show help for unknown subcommands
        parser.print_help()
