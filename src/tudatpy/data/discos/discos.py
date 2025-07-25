import requests

class DiscosQuery:
    def __init__(self, token, url="https://discosweb.esoc.esa.int"):
        """
        Initializes the DiscosQuery object with the provided token and base URL.
        """
        self.token = token
        self.url = url
        self.api_version = '2'
        self.headers = {
            'Authorization': f'Bearer {self.token}',
            'DiscosWeb-Api-Version': self.api_version,
        }

    def query_object(self, norad_id, verbose=True):
        """
        Queries the DISCOS database for an object with the specified NORAD ID.

        Parameters:
        - norad_id (int or str): NORAD catalog ID of the satellite
        - verbose (bool): whether to print the results (default True)

        Returns:
        - dict: attributes of the queried object (e.g. mass, launch date, etc.)
        """
        response = requests.get(
            f'{self.url}/api/objects?filter=eq(satno,{norad_id})',
            headers=self.headers,
        )

        if response.ok:
            attributes = response.json()['data'][0]['attributes']
            if verbose:
                print(attributes)
            return attributes
        else:
            errors = response.json().get('errors', 'Unknown error')
            if verbose:
                print(errors)
            return None

    def get_object_attribute(self, norad_id, attribute, verbose=True):
        attributes = self.query_object(norad_id)
        if attributes:
            if attribute in attributes:
                return attributes[attribute]
            else:
                print('Attribute not found. Exiting...\n')
                exit()


