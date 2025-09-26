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
            f'{self.url}/api/objects/{norad_id}',
            headers=self.headers
        )

        if response.ok:
            attributes = response.json()['data']['attributes']
            if verbose:
                print(attributes)
                print(attributes.get('mass', 'Mass not available'))
            return attributes
        else:
            errors = response.json().get('errors', 'Unknown error')
            if verbose:
                print(errors)
            return None