import os
import sys
import tempfile
import shutil
from datetime import date


HELP_TEXT = \
"""One of the following arguments must be supplied:
- 'major': to increase major version number, e.g. from 1.2.3 to 2.0.0
- 'minor': to increase minor version number, e.g. from 1.2.3 to 1.3.0
- 'patch': to increase patch number, e..g from 1.2.3 to 1.2.4"""

LEVELS = ['major', 'minor', 'patch']


def main(argv):
    try:
        level_index = LEVELS.index(argv[0].lower())  # 0, 1 or 2
    except (IndexError, ValueError):
        print(HELP_TEXT)
        sys.exit()

    old_version = existing_version()

    new_version = [0, 0, 0]
    for index, value in enumerate(new_version[0:level_index]):
        new_version[index] = old_version[index]
    new_version[level_index] = old_version[level_index] + 1

    print("Previous version: {}.{}.{}".format(*old_version))
    print("Current version:  {}.{}.{}".format(*new_version))

    update_package_setup(new_version)
    print("Packup setup.py updated.")

    update_doc_conf(new_version)
    print("Documentation conf.py updated.")

    update_changelog(new_version)
    print("CHANGELOG.txt updated.")


def existing_version():
    with open('setup.py') as setup_file:
        for line in setup_file:
            if line.strip().startswith('version'):
                return [int(s) for s in line.strip().split('=')[1].strip("',").split('.')]


def replace_file_content(file_name, content):
    with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', delete=False) as temp_file:
        for line in content:
            temp_file.write(line)
    shutil.copy(temp_file.name, file_name)
    os.remove(temp_file.name)


def update_package_setup(new_version):
    file_name= 'setup.py'

    new_content = []
    with open(file_name) as f:
        for line in f:
            if line.strip().startswith('version'):
                line = "    version='{}.{}.{}',\n".format(*new_version)
            new_content.append(line)

    replace_file_content(file_name, new_content)


def update_doc_conf(new_version):
    file_name = './docs/source/conf.py'

    new_content = []
    with open(file_name) as f:
        for line in f:
            if line.strip().startswith('version'):
                line = "version = '{}.{}'\n".format(*new_version[0:2])
            elif line.strip().startswith('release'):
                line = "release = '{}.{}.{}'\n".format(*new_version)
            new_content.append(line)

    replace_file_content(file_name, new_content)


def update_changelog(new_version):
    file_name = 'CHANGELOG.txt'

    header = 'version {}.{}.{} ({})'.format(new_version[0], new_version[1], new_version[2], date.today())
    new_content = [header + '\n', '-' * len(header) + '\n']
    with open(file_name) as f:
        for line in f:
            new_content.append(line)

    replace_file_content(file_name, new_content)


if __name__ == "__main__":
    main(sys.argv[1:])