import sys


def job(smiles):
    # Task: Write database code here using smiles
    pass


if __name__ == "__main__":
    if len(sys.argv) == 2:
        # For Directory

        # Task: Write Smiles retrieving code from recently uploaded sdf file from upload folder
        smiles = ''
        job(smiles)
    elif len(sys.argv) == 3:
        # For Specific file

        # Task: Write Smiles retrieving code from recently uploaded sdf file from upload specific file
        smiles = ''
        job(smiles)
    else:
        pass
        # Logging code
