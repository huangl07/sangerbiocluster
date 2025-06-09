# -*- coding:utf-8 -*-
"""
Last-edit: 2023/5/9
Author: yiwei.tang
mail: yiwei.tang@majorbio.com
"""
import argparse
import shutil
import os
from pathlib import Path
import re
import yaml


def publish_files(base_dir: str, pub_dir: str, config_data: list,
                  pid: int) -> list:
    """
    按需整理文件，基于publish_dir关键词from/to
    """
    base = Path(base_dir)
    pub = Path(pub_dir)
    info_list = []
    for index, config in enumerate(config_data):
        check = True
        if "check" in config.keys():
            check = base.joinpath(config["check"]).exists()
        if not check:
            print("rule {} check failed!".format(index))
            continue
        if "mkdir" in config.keys():
            mkd = config["mkdir"].format(id=str(pid).zfill(2))
            pub.joinpath(mkd).mkdir()
            print(config["mkdir"] + " created!")
            info_list.append([
                str(Path(mkd)),
                Path(mkd).name, config["desc"],
                "folder" if pub.joinpath(mkd).parent != pub else "basefolder"
            ])
            continue
        if "from" not in config.keys():
            config["from"] = ""
        if "to" not in config.keys():
            continue
        if not base.joinpath(config["from"]).exists():
            print(config["from"] + " does not exist!")
            continue
        if "pattern" not in config.keys():
            config["pattern"] = ".*"
        from_list = [
            i.name for i in base.joinpath(config["from"]).glob("*")
            if re.compile(config["pattern"]).match(i.name)
        ]
        if len(from_list) == 0:
            print(config["from"] + '/' + config["pattern"] + ' not found!')
        from_list.sort()
        for src_file in from_list:
            match_obj = re.compile(config["pattern"]).match(src_file)
            if not match_obj:
                continue
            to_path = config["to"].format(id=str(pid).zfill(2),
                                          **match_obj.groupdict())
            if not pub.joinpath(to_path).exists():
                print(config["to"] + " has not been created!")
                break
            os.system("cp -rl " +
                      str(base.joinpath(config["from"]).joinpath(src_file)) +
                      " " + str(pub.joinpath(to_path)))
            info_list.append([
                str(Path(to_path) / src_file), src_file,
                config["desc"].format(id=str(pid).zfill(2),
                                      **match_obj.groupdict()), "file"
            ])
            print(
                str(base.joinpath(config["from"]).joinpath(src_file)) +
                " => " + to_path)
    return info_list


def read_yaml(path: str) -> dict:
    """
    读取yaml配置文件
    """
    with open(path, 'r', encoding="utf-8") as f_in:
        content = f_in.read()
    config = yaml.load(content, Loader=yaml.Loader)
    return config["data"]


def generate_body_html(path: str, body_info: list):
    """
    生成body.html
    """
    with open(path, "w", encoding="utf-8") as body_out:
        for info in body_info:
            if info[3] == "basefolder":
                body_out.write('''<tr data-tt-id='{filepath}'>
                    <td><span class='folder'>{filename}</span></td>
                    <td>{filedesc}</td>
                    <td>文件夹</td>
                    </tr>\n\n'''.format(filepath=info[0].replace("/", "-"),
                                        filename=info[1],
                                        filedesc=info[2]))
            elif info[3] == "folder":
                parentfile = str(Path(info[0]).parent).replace("/", "-")
                body_out.write(
                    '''<tr data-tt-id='{filepath}' data-tt-parent-id='{parentfile}'>
                    <td><span class='folder'>{filename}</span></td>
                    <td>{filedesc}</td>
                    <td>文件夹</td>
                    </tr>\n\n'''.format(filepath=info[0].replace("/", "-"),
                                        filename=info[1],
                                        filedesc=info[2],
                                        parentfile=parentfile))
            else:
                parentfile = str(Path(info[0]).parent).replace("/", "-")
                body_out.write(
                    '''<tr data-tt-id='{fileid}' data-tt-parent-id='{parentfile}'>
                    <td><span class='file'><a href='{filepath}'>{filename}</a></span></td>
                    <td>{filedesc}</td>
                    <td>文件</td>
                    </tr>\n\n'''.format(fileid=info[0].replace("/", "-"),
                                        filepath=info[0],
                                        filename=info[1],
                                        filedesc=info[2],
                                        parentfile=parentfile))


if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="根据yaml配置文件重排输出文件")
    parse.add_argument("-c", "--config", help="yaml配置文件", required=True)
    parse.add_argument("-i", "--indir", help="输入文件夹", required=True)
    parse.add_argument("-o", "--outdir", help="输出文件夹", required=True)
    parse.add_argument("-t",
                       "--clean",
                       action="store_true",
                       help="是否清空输出文件夹",
                       default=False)
    args = parse.parse_args()
    data = read_yaml(args.config)
    print("Base Dir: " + args.indir)
    print("Publish Dir: " + args.outdir)
    if args.clean:
        shutil.rmtree(args.outdir)
    if not Path(args.outdir).exists():
        Path(args.outdir).mkdir()
    body_list = []
    PID = 0
    for basefolder, data_pieces in data.items():
        if Path(args.indir).joinpath(basefolder).exists():
            PID += 1
            body_list += publish_files(args.indir, args.outdir, data_pieces,
                                       PID)
    generate_body_html("body.html", body_list)
