# https://gist.github.com/alanbchristie/ce49248560060c9562d7f7e34f8172af

# I have never got this to work.

import os
import base64
import requests
import functools
from typing import Optional
class FragalysisCompoundSetSubmitter:
    def __init__(self,
                 fragalysis_hostname: str="fragalysis.diamond.ac.uk",
                 method: str = "viewer/upload_cset",
                 settings: Optional[dict] = None):
        self.url = f"https://{fragalysis_hostname}/{method}/"
        self.settings = {} if settings is None else settings
        self.keycloak_url = self.get('keycloak_url')
        self.keycloak_realm = self.get('keycloak_realm', 'xchem')
        self.keycloak_client_id = self.get('keycloak_client_id')
        self.keycloak_client_secret = self.get('keycloak_client_secret')
        self.keycloak_username = self.get('keycloak_username')
        self.keycloak_password = self.get('keycloak_password')



    def submit_cset(self, target_name, sdf_path, update_set='None', dry_run=False,
                    pdb_zip_path=None,
                    add=False) -> str:
        """

        :param target_name:
        :param sdf_path:
        :param update_set: "".join(submitter_name.split()) + '-' + "".join(method.split())
        :param dry_run:
        :param pdb_zip_path:
        :param add:
        :return:
        """
        if not add:
            payload = {'target_name': target_name,
                       'submit_choice': 0 if dry_run else 1,
                       'update_set': update_set}
        else:
            payload = {'target_name': target_name,
                       'submit_choice': 0 if dry_run else 1,
                       'update_set': update_set}

        response: requests.Response = requests.post(self.url, headers=self.headers, data=payload)
        response.raise_for_status()
        assert response is not None
        return response.text

    def get(self, key: str, default: str = None):
        if key in self.settings:
            return self.settings[key]
        elif key.upper() in self.settings:
            return self.settings[key.upper()]
        elif key in os.environ:
            return os.environ[key]
        elif key.upper() in os.environ:
            return os.environ[key.upper()]
        elif default is not None:
            return default
        else:
            raise ValueError(f"Environment variable '{key}' not found")

    @functools.cached_property
    def keycloak_access_token(self) -> str:
        """Gets an 'access token' from Keycloak.
        If successful we'll get the token (a big long string)
        """
        realm_url: str = f"{self.keycloak_url}/realms/{self.keycloak_realm}"
        url = f"{realm_url}/protocol/openid-connect/token"
        data = (
            f"client_id={self.keycloak_client_id}"
            f"&grant_type=password"
            f"&username={self.keycloak_username}"
            f"&password={self.keycloak_password}"
            f"&client_secret={self.keycloak_client_secret}"
        )
        headers = {"Content-Type": "application/x-www-form-urlencoded"}
        response : requests.Response = requests.post(url, headers=headers, data=data, timeout=4.0)
        response.raise_for_status()
        assert "access_token" in response.json()
        return response.json()["access_token"]

    @functools.cached_property
    def csrf_token(self) -> str:
        client = requests.session()
        # Retrieve the CSRF token first
        client.get(self.url)  # sets cookie
        if 'csrftoken' in client.cookies:
            # Django 1.6 and up
            csrftoken = client.cookies['csrftoken']
        else:
            # older versions
            csrftoken = client.cookies['csrf']
        return csrftoken

    @functools.cached_property
    def headers(self) -> dict:
        if self.keycloak_client_secret:
            return {"X-CSRFToken": self.csrf_token,
                       "Cookie": f"csrftoken={self.csrf_token}",
                       "Authorization": f"Bearer {self.keycloak_access_token}"}
        else:
            userpass = f"{self.keycloak_username}:{self.keycloak_password}"
            basic = str(base64.b64encode(userpass.encode("utf-8")), "utf-8")
            return {"X-CSRFToken": self.csrf_token,
                       "Cookie": f"csrftoken={self.csrf_token}",
                       "Authorization": f"Basic {basic}"}