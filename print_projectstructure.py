# %%
import os

# a. Plain style
def print_lines(startpath, output_file, ignore_dirs=None, ignore_files=None):
    if ignore_dirs is None:
        ignore_dirs = []
    if ignore_files is None:
        ignore_files = []

    with open(output_file, 'w', encoding='utf-8') as f:
        for root, dirs, files in os.walk(startpath):
            # Filter out ignored directories and files
            dirs[:] = [d for d in dirs if d not in ignore_dirs]
            files[:] = [f for f in files if f not in ignore_files]
            # Sort for consistent output
            dirs.sort()
            files.sort()

            level = root.replace(startpath, '').count(os.sep)
            indent = '    ' * level
            print(f"{indent}{os.path.basename(root)}/")
            subindent = '    ' * (level + 1)
            for f in files:
                print(f"{subindent}{f}")

# b. A little prettier style
def print_tree(startpath, output_file, ignore_dirs=None, ignore_files=None):
    if ignore_dirs is None:
        ignore_dirs = []
    if ignore_files is None:
        ignore_files = []
    
    with open(output_file, 'w', encoding='utf-8') as f:
        for root, dirs, files in os.walk(startpath):
            dirs[:] = [d for d in dirs if d not in ignore_dirs]
            files[:] = [f for f in files if f not in ignore_files]
            dirs.sort()
            files.sort()

            level = root.replace(startpath, '').count(os.sep)
            indent = '│   ' * level
            branch = '├── ' if level > 0 else ''
            f.write(f"{indent}{branch}{os.path.basename(root)}/\n")

            for i, file in enumerate(files):
                subindent = '│   ' * (level + 1)
                connector = '└── ' if i == len(files) - 1 else '├── '
                f.write(f"{subindent}{connector}{file}\n")

# Example usage
project_root = os.path.abspath('.')
output_txt = os.path.join(project_root, 'project_structure.txt')
folders_toignore = ['.git', '.Rproj.user', '.venv']
files_toignore = ['.gitignore', 'git_project.Rproj']
print_tree(project_root, output_txt, ignore_dirs=folders_toignore, ignore_files=files_toignore)
# End