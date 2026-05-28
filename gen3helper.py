"""
Helper client for interfacing with Gen3 APIs
"""

from gen3.auth import Gen3Auth
from gen3.file import Gen3File
from gen3.tools.download.drs_download import DownloadManager, Downloadable
from gen3.submission import Gen3Submission
import os


class gen3Client:
    def __init__(self, endpoint: str, credential_file: str = None):
        self.endpoint = endpoint
        self.hostname = endpoint.removeprefix("https://")
        self.auth = Gen3Auth(endpoint=self.endpoint, refresh_file=credential_file)
        self.file_client = Gen3File(self.auth)
        self.submission_client = Gen3Submission(self.auth)

    def download_files(
        self,
        file_ids: list[str],
        save_directory="./Downloads",
        show_progress: bool = True,
    ):
        """Download a list of files from Gen3"""
        download_list = [Downloadable(object_id=f) for f in file_ids]
        manager = DownloadManager(
            hostname=self.hostname,
            auth=self.auth,
            download_list=download_list,
            show_progress=show_progress,
        )
        os.makedirs(save_directory, exist_ok=True)
        existing_files = os.listdir(save_directory)
        manager.download_list = list(
            filter(lambda x: x.file_name not in existing_files, manager.download_list)
        )
        return manager.download(
            manager.download_list,
            save_directory=save_directory,
            show_progress=show_progress,
        )
    def get_project_files(self, project_id):
        """Get all files from a project"""
      
        # Common file node types to try
        file_node_types = [
            'submitted_aligned_reads',
            'aligned_reads', 
            'submitted_unaligned_reads',
            'simple_germline_variation',
            'reference_file',
            'core_metadata_collection',
            'submitted_somatic_mutation'
        ]
        
        all_files = []
        
        for node_type in file_node_types:
            query = f"""
            {{
                {node_type}(project_id: "{project_id}", first: 1000) {{
                    id
                    object_id
                    file_name
                    file_size
                    data_format
                    data_type
                    md5sum
                    project_id
                }}
            }}
            """
            
            try:
                response = self.submission_client.query(query)
                files = response.get('data', {}).get(node_type, [])
                if files:
                    print(f"Found {len(files)} files in {node_type}")
                    for f in files:
                        f['node_type'] = node_type
                    all_files.extend(files)
            except Exception as e:
                # Skip node types that don't exist
                continue
        
        return all_files

    def list_open_projects(self):
        query = """
        {
            project(first: 100) {
                project_id
                name
                state
                code
            }
        }
        """
        response = self.submission_client.query(query)
        print("Open access projects:")
        print(response['data']['project'])
        print("-"*100)

