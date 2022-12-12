import argparse
import os
import os.path
import tempfile

import tqdm

from apiclient.http import MediaIoBaseDownload
from google.oauth2 import service_account
from googleapiclient.discovery import build

TOP_DIR_NAME = 'IRIS_data'
CHUNK_SIZE = 1024 * 1024 * 8  # 8 MB
SCOPES = ['https://www.googleapis.com/auth/drive']


def parse_args():
    parser = argparse.ArgumentParser(
        description=('download IRIS files from google drive'))
    parser.add_argument('--iris-folder-id',
                        help=('ID of IRIS_data folder on google drive'
                              ' (can be found in download url)'))
    parser.add_argument(
        '--dest-dir',
        help='path to directory where IRIS_data/ will be written')
    parser.add_argument('--download-all',
                        action='store_true',
                        help='download all data')
    parser.add_argument('--list-files',
                        action='store_true',
                        help='write available files to --selected-tsv')
    parser.add_argument('--selected-tsv',
                        help='path to a tsv with selected files to download')
    parser.add_argument(
        '--api-key-json-path',
        required=True,
        help='path to the .json file that is the service account key')

    args = parser.parse_args()
    if args.download_all:
        if not (args.iris_folder_id and args.dest_dir):
            parser.error(
                '--download-all requires --iris-folder-id and --dest-dir')
    elif args.list_files:
        if not (args.iris_folder_id and args.selected_tsv):
            parser.error(
                '--list-files requires --iris-folder-id and --selected-tsv')
    elif not (args.dest_dir and args.selected_tsv):
        parser.error(
            'download specific files with --dest-dir and --selected-tsv.'
            ' Otherwise use --download-all or --list-files')

    return args


def main():
    args = parse_args()

    credentials = service_account.Credentials.from_service_account_file(
        args.api_key_json_path, scopes=SCOPES)
    with build('drive', 'v3', credentials=credentials) as drive_service:
        if args.download_all:
            download_all_files(args.iris_folder_id, args.dest_dir,
                               drive_service)
            return

        if args.list_files:
            list_all_files(args.iris_folder_id, args.selected_tsv,
                           drive_service)
            return

        download_selected_files(args.dest_dir, args.selected_tsv,
                                drive_service)


def list_files_recursive(parent_id, drive_service):
    results = list()
    files_c = drive_service.files()
    request = files_c.list(q="'{}' in parents".format(parent_id))

    response = request.execute()
    found = response.get('files', list())
    for file_dict in found:
        file_id = file_dict['id']
        file_name = file_dict['name']
        is_folder = 'folder' in file_dict['mimeType']
        if is_folder:
            sub_results = list_files_recursive(file_id, drive_service)
            results.append({
                'folder': file_name,
                'id': file_id,
                'files': sub_results
            })
        else:
            results.append({'name': file_name, 'id': file_id})

    return results


def write_tsv_line(columns, tsv_handle):
    tsv_handle.write('{}\n'.format('\t'.join(columns)))


def write_file_tsv(all_files, parent_dir_path, tsv_handle):
    files = list()
    folders = list()
    for file_dict in all_files:
        if 'folder' in file_dict:
            folders.append(file_dict)
        else:
            files.append(file_dict)

    files.sort(key=lambda d: d['name'])
    folders.sort(key=lambda d: d['folder'])
    for file_dict in files:
        full_path = os.path.join(parent_dir_path, file_dict['name'])
        write_tsv_line([full_path, file_dict['id']], tsv_handle)

    for folder_dict in folders:
        full_path = os.path.join(parent_dir_path, folder_dict['folder'])
        write_file_tsv(folder_dict['files'], full_path, tsv_handle)


def download_all_files(iris_folder_id, dest_dir, drive_service):
    all_files = list_files_recursive(iris_folder_id, drive_service)

    temp_name = None
    try:
        with tempfile.NamedTemporaryFile(delete=False) as temp_handle:
            temp_name = temp_handle.name

            write_file_tsv(all_files, TOP_DIR_NAME, temp_handle)

        download_selected_files(dest_dir, temp_name, drive_service)
    finally:
        if temp_name:
            os.remove(temp_name)


def list_all_files(iris_folder_id, selected_tsv, drive_service):
    all_files = list_files_recursive(iris_folder_id, drive_service)
    with open(selected_tsv, 'wt') as tsv_handle:
        write_file_tsv(all_files, TOP_DIR_NAME, tsv_handle)


def download_file(file_id, dest_path, drive_service):
    files_c = drive_service.files()
    request = files_c.get_media(fileId=file_id)
    dir_path = os.path.dirname(dest_path)
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)

    with open(dest_path, 'wb') as out_handle:
        downloader = MediaIoBaseDownload(out_handle,
                                         request,
                                         chunksize=CHUNK_SIZE)
        progress = tqdm.tqdm(desc=dest_path, total=1.0, unit='file')
        progress_so_far = 0
        done = False
        while done is False:
            status, done = downloader.next_chunk()
            if status:
                new_progress = status.progress()
                additional = new_progress - progress_so_far
                progress_so_far = new_progress
                progress.update(additional)

        progress.close()


def download_selected_files(dest_dir, selected_tsv, drive_service):
    with open(selected_tsv, 'rt') as tsv_handle:
        for line in tsv_handle:
            columns = line.strip().split('\t')
            if len(columns) != 2:
                raise Exception(
                    'expected 2 columns in {}'.format(selected_tsv))

            drive_file_path = columns[0]
            file_id = columns[1]

            local_file_path = os.path.join(dest_dir, drive_file_path)
            download_file(file_id, local_file_path, drive_service)


if __name__ == '__main__':
    main()
