from pathlib import Path
import os
import shutil


def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


if __name__ == "__main__":

    # Delete all headers in /src
    for path in Path('tests').rglob('*.dat'):
        print(path.name)
        if path.name != "CMakeLists.txt":
            print(path)
        # os.remove(str(path))
        # pass
        # print(str(path))
            path_rel_test_dest = os.path.join(*str(path).split("/")[:])
            dest = os.path.join("tests", "data", path.name)
            print(path, " to ", dest)
            try:
                shutil.copyfile(path, dest)
            except shutil.SameFileError:
                pass

            # print(dest)
    #     try:
    #         shutil.copyfile(path, dest)
    #     except FileNotFoundError:
    #         pass
    #     # print(path_rel_test_dest)
    # #     if os.path.isdir(path):
    # #         print("-----" * 10)
    # #         print(path)
    #         print(os.listdir(path))
    # #
    # # Extract all tests
    # for path in Path('tests').rglob('*.h'):
    #         # if os.path.isdir(path):
    #
    #         path_rel_test_dest = os.path.join(*str(path).split("/")[1:-1])
    #         print(path_rel_test_dest)
    #         dest = os.path.join('tests', "include", "tudat", path_rel_test_dest)
    #         src = path
    #
    #         print(src, ' to ', dest)
    #         if os.path.exists(dest):
    #             pass
    #         else:
    #             os.makedirs(dest)
    #         shutil.copyfile(src, os.path.join(dest, path.name))
    #
    #         # os.makedirs(dest)
    #         # print(path, ' to ', dest)
    #         # copytree(src, dest)
