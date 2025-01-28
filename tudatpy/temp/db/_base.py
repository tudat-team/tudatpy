import urllib.request
import os
from clint.textui import progress
import requests
import json


class ResourceTool:

    def __init__(self, resource_url, resource_dir, resource_name, extension):
        self._resource_url = resource_url
        self._resource_dir = resource_dir
        self._resource_name = resource_name
        self._extension = extension
        self._encoding = "utf-8"
        self._html = None

    def _get_destined_path(self):
        return os.path.join(self._resource_dir,
                            self._resource_name + self._extension)

    def _get_file(self):
        if not self.is_downloaded:
            self.download()
        return open(self._get_destined_path(), "r")

    def download(self):
        print(f"Downloading {os.path.basename(self._get_destined_path())}")
        dirname = os.path.dirname(self._get_destined_path())
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        r = requests.get(self._resource_url, stream=True)
        with open(self._get_destined_path(), 'wb') as f:
            total_length = int(r.headers.get('content-length'))
            for chunk in progress.bar(
                    r.iter_content(chunk_size=1024),
                    expected_size=(total_length / 1024) + 1):
                if chunk:
                    f.write(chunk)
                    f.flush()

    @property
    def is_downloaded(self):
        return os.path.exists(self._get_destined_path())

    @property
    def string(self):
        return self._get_file().read()

    @property
    def path(self, ensure_downloaded=True):
        if not self.is_downloaded and ensure_downloaded:
            self.download()
        return self._get_destined_path()


class DataBaseIndex:

    def __init__(self, dict, resource_root="./cache"):
        self._name = dict["name"]
        self._resources = dict["resources"]
        self._resource_root = resource_root

    def __getattr__(self, key):
        if self._resources[key]["type"] == "file":
            return _BaseResource(
                resource_url=self._resources[key]["resource_url"],
                resource_dir=os.path.join(self._resource_root, self._name),
                resource_name=os.path.basename(self._resources[key]["resource_url"]).split(".")[0],
                extension="." + os.path.basename(self._resources[key]["resource_url"]).split(".")[1]
            )
        elif self._resources[key]["type"] == "directory":
            _key_list = []
            _cls_list = []
            for item in self._resources[key]["items"]:
                file_name = os.path.basename(item["resource_url"])
                _key_list.append(file_name)
                _cls_list.append(_BaseResource(
                    resource_url=item["resource_url"],
                    resource_dir=os.path.join(self._resource_root, self._name, key),
                    resource_name=file_name.split(".")[0],
                    extension="." + file_name.split(".")[1]
                ))
            return dict(zip(_key_list, _cls_list))
        else:
            raise KeyError("type must be specified for a resource.")

    @classmethod
    def from_json(cls, path, resource_root="~/.cache"):
        with open(path, 'r') as j:
            contents = json.loads(j.read())
        return cls(contents, resource_root)

    @classmethod
    def from_dict(cls, dict, resource_root="~/.cache"):
        return cls(dict, resource_root)

    def to_json(self, path, root_dir="~/.cache"):



if __name__ == "__main__":

    celestrak_db = DataBaseIndex.from_json("celestrak.json")
    print(celestrak_db.active_satellites.path)
    print(celestrak_db.space_stations.path)

    tudat_space_db = DataBaseIndex.from_json("tudat_space.json")
    print(tudat_space_db.spice_kernels["de430_mar097_small.bsp"].path)
    print(tudat_space_db.spice_kernels["naif0012.tls"].path)
    print(tudat_space_db.spice_kernels["tudat_merged_spk_kernel.inp"].path)
    print(tudat_space_db.spice_kernels["moon_assoc_pa.tf"].path)


