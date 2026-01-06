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

    def query_object(self, sat_id, is_discos_id=False, verbose=True) -> dict[str] | None:
        """
        Queries the DISCOS database using either a NORAD ID (default) or DISCOS ID.

        Parameters
        ----------
        sat_id: str
            The ID of the satellite to query.
        is_discos_id: bool
            If True, treats sat_id as an internal DISCOS ID.
            If False (default), treats sat_id as a NORAD ID (satno).

        verbose: bool
            Whether to print the results (default True).

        Returns
        ----------
        dict[str] | None
            The attributes of the queried object.
        """

        if is_discos_id:
            query_url = f'{self.url}/api/objects/{sat_id}'
        else:
            query_url = f'{self.url}/api/objects?filter=eq(satno,{sat_id})'

        response = requests.get(query_url, headers=self.headers)

        if response.ok:
            data = response.json().get('data')

            # Handle API response structure:
            if isinstance(data, list):
                if not data:
                    if verbose: print(f"No object found for queried satellite.")
                    return None
                attributes = data[0]['attributes']
            else:
                attributes = data['attributes']

            if verbose:
                print(attributes)

            return attributes

        else:
            # Handle API Errors
            errors = response.json().get('errors', 'Unknown error')
            if verbose:
                print(f"API Error: {errors}")
            return None