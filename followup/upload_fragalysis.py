"""

Copied from https://github.com/xchem/fragalysis-api/blob/b09a4b978c2e2bdfdc2d9d0b4c85178b24e949c7/fragalysis_api/xcextracter/computed_set_update.py


EXAMPLES:
=========
---- NB: The major difference here is the REQ_URL. For new data, or to overwrite data, use https://fragalysis.diamond.ac.uk/viewer/upload_cset/.
         To add new molecules to an existing set, use https://fragalysis.diamond.ac.uk/viewer/update_cset/. ----

to overwrite an existing cset:
------------------------------
taskurl = update_cset(REQ_URL='https://fragalysis.diamond.ac.uk/viewer/upload_cset/',
                      target_name='Mpro',
                      submit_choice='1',
                      upload_key='1',
                      update_set='WT-xCOS3-ThreeHop',
                      sdf_path='/Users/uzw12877/Downloads/Test_upload/Top_100_three_hop_XCOS_1.4_2020-07-28.sdf',
                      pdb_zip_path='/Users/uzw12877/Downloads/Test_upload/receptor_pdbs.zip')
task_response, json_results = get_task_response(taskurl)
print(task_response)
print(json_results)

to update an existing cset:
---------------------------
taskurl = update_cset(REQ_URL='https://fragalysis.diamond.ac.uk/viewer/update_cset/',
                      target_name='Mpro',
                      update_set='WT-xCOS3-ThreeHop',
                      sdf_path='/Users/uzw12877/Downloads/Test_upload/Top_100_three_hop_XCOS_1.4_2020-07-28 copy.sdf',
                      pdb_zip_path='/Users/uzw12877/Downloads/Test_upload/receptor_pdbs copy.zip',
                      add=True)
task_response, json_results = get_task_response(taskurl)
print(task_response)
print(json_results)

"""

import requests
import time
import logging

logger = logging.getLogger(__name__)

REQ_URL = 'https://fragalysis.diamond.ac.uk/viewer/upload_cset/'


def get_csrf(REQ_URL=REQ_URL):
    """Get a csrf token from the request url to authenticate further requests
    Parameters
    ----------
    REQ_URL string
        The URL that you want to make a request against after getting the token

    Returns
    -------
    csrftoken
        csrf token to use for further requests against the same URL
    """
    client = requests.session()
    # Retrieve the CSRF token first
    client.get(REQ_URL)  # sets cookie
    if 'csrftoken' in client.cookies:
        # Django 1.6 and up
        csrftoken = client.cookies['csrftoken']
    else:
        # older versions
        csrftoken = client.cookies['csrf']
    return csrftoken


def update_cset(target_name, sdf_path, update_set='None', submit_choice=None, upload_key=None,
                pdb_zip_path=None, add=False, REQ_URL=REQ_URL):
    """Send data to <root_url>/viewer/upload_cset/ to overwrite an existing computed set, or to
    <root_url>/viewer/update_cset/ to add new molecules without deleting the old ones.

    Parameters
    ----------
    REQ_URL: str
        request URL for the upload (e.g. https://fagalysis.diamond.ac.uk/viewer/upload_cset/ or viewer/update_cset)
    target_name: str
        the name of the target in Fragalysis that the computed set is for
    update_set: str
        the name of the computed set you want to update
        (can be found with: "".join(submitter_name.split()) + '-' + "".join(method.split()),
        where submitter_name is the name in the submitter_name field in the blank mol of the uploaded sdf file,
        and method is the method field in the blank mol of the uploaded sdf file). Leave blank if you are adding
        a set for the first time
    sdf_path: str
        path to the sdf file to upload
    submit_choice: int
        0 for validate, 1 for upload (not required for update - ie. viewer/update_cset)
    upload_key: str
        upload key, not currently turned on, so can be any value, but not blank or null (optional)
    pdb_zip_path: str
        path to the zip file of pdb's to upload (optional)
    add: bool
        set to True if updating a computed set without overwriting it completely (for <root_url>/viewer/update_cset/)

    Returns
    -------
    taskurl str
        the URL to check for the status of the upload
    """
    logger.debug(f'Submitting files to update {update_set}...')
    csrf_token = get_csrf(REQ_URL)
    if not add:
        payload = {'target_name': target_name,
                   'submit_choice': submit_choice,
                   'upload_key': upload_key,
                   'update_set': update_set}
    else:
        payload = {'target_name': target_name,
                   'update_set': update_set}
    files = [('sdf_file', open(sdf_path, 'rb')),]
    if pdb_zip_path:
        files.append(('pdb_zip', open(pdb_zip_path, 'rb')))
    headers = {'X-CSRFToken': csrf_token,
               'Cookie': f'csrftoken={csrf_token}'}
    response = requests.post(REQ_URL, headers=headers, data=payload, files=files, timeout=600)
    lines = response.text.split('\n')
    taskurl = None
    for l in lines:
        if 'taskUrl = "/viewer/upload_task/' in l:
            taskid = l.split('/')[-2]
            logger.debug(f'upload task id: {taskid}')
            taskurl = f'{REQ_URL.replace("/viewer/upload_cset/", "/viewer/upload_task/")}{taskid}'
        elif 'taskUrl = "/viewer/update_task/' in l:
            taskid = l.split('/')[-2]
            logger.debug(f'upload task id: {taskid}')
            taskurl = f'{REQ_URL.replace("/viewer/update_cset/", "/viewer/update_task/")}{taskid}'
            if taskurl:
                break
    assert taskurl, f'Something went wrong with the upload/update request! \
                        Please try again or email frank.vondelft@cmd.ox.ac.uk for help.\
                        Response: {response.text}'
    logger.debug('pinging task to check status...')
    requests.get(taskurl)
    delta_t = 2
    for trial in range(600 // delta_t):
        task_response: requests.Response = requests.get(taskurl, timeout=600)
        data = task_response.json()
        status = data.get('upload_task_status', '') + data.get('update_task_status', '')
        if status == "SUCCESS" or status == "FAILURE":
            return status, task_response.json()
        time.sleep(delta_t)
