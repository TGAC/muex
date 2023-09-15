import defopt
from pathlib import Path

# TODO allow relative paths?
# TODO print mappings to file

def main(indir: str, outdir: str):
    """
    Creates consistently-named symlinks in `outdir` to all files in `indir`. Files are named `file0.fq.gz`,`file1.fq.gz`... in alphabetical order of their names in `indir`.

    :param indir: Absolute path to directory containing input files
    :param outdir: Absolute path to directory to contain desired symlinks
    """
    inpath = Path(indir)
    if "~" in indir:
        inpath = inpath.expanduser()
    
    outpath = Path(outdir)
    if "~" in indir:
        outpath = outpath.expanduser()

    outpath.mkdir(parents=True, exist_ok=True)

    # Path.resolve() to get absolute path (replaces symlinks)
    files = sorted([x.resolve() for x in list(inpath.iterdir()) if x.is_file()])
    num_files = len(files)

    for i in range(num_files):
        symlink_filepath = outpath / f"file{i}.fq.gz"
        symlink_filepath.symlink_to(files[i])

    return

if __name__ == "__main__":
    defopt.run(main)