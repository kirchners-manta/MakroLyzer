import os


def output_path(args, filename):
    """Place relative output filenames under the requested directory."""
    directory = args.get('output_dir')
    if directory is None:
        return filename
    return os.path.join(directory, filename)
