from gen3helper import gen3Client

COMMONS = 'https://gen3.biodatacatalyst.nhlbi.nih.gov'
CRED_FILE = '/Users/buzgalbraith/.gen3/credentials.json'
target_project = 'tutorial-synthetic_data_set_1'


if __name__ == "__main__":
    client = gen3Client(endpoint=COMMONS, credential_file=CRED_FILE)
    client.list_open_projects()
    project_files = client.get_project_files(target_project)
    oid = 'dg.4503/0fbb8b5d-81a5-4928-a42d-7cac707f746e'
    client.download_files(file_ids=[oid])
    [x for x in project_files if 'ALL.chr17' in x['file_name']]