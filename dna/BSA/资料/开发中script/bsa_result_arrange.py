# -*- coding:utf-8 -*-
"""
LastEditTime: 2023/06/06
Author: yuan.xu
mail: yuan.xu@majorbio.com
"""

import argparse
import shutil
import os
import glob
from pathlib import Path
import re
import yaml
from bin.bsa_function import arrange_pop_table



def read_yaml(path: str) -> dict:
    """
    读取yaml配置文件
    """
    with open(path, 'r') as f_in:
        content = f_in.read()
    config = yaml.load(content, Loader=yaml.Loader)
    return config["arrange"]


def publish_files(base_dir: str, pub_dir: str, config_data: list) -> list:
    """
    按需整理文件,基于publish_dir关键词from/to
    """
    base = Path(base_dir)
    pub = Path(pub_dir)
    for index, config in enumerate(config_data):
        check = True
        if "check" in config.keys():
            check = base.joinpath(config["check"]).exists()
        if not check:
            print("rule {} check failed!".format(index))
            continue
        if "mkdir" in config.keys():
            if pub.joinpath(config["mkdir"]).exists():
                continue
            else:
                pub.joinpath(config["mkdir"]).mkdir()
                print(config["mkdir"] + " created!")
                continue
        if "file" in config.keys():
            if isinstance(config["file"], str):
                from_list = [
                    i.name for i in base.joinpath(config["check"]).glob("*")
                    if re.compile(config["file"]).match(i.name)
                    ]
                if len(from_list) == 0:
                    print("{} not found!".format(config["file"]))
                elif len(from_list) == 1:
                    source_file = str(base.joinpath(config["check"]).joinpath(from_list[0]))
                    dest_file = str(pub.joinpath(config["to"]).joinpath(config["new_file"]))
                    shutil.copyfile(source_file, dest_file)
                    print("文件 {} 拷贝成 {} 成功".format(source_file, dest_file))
                else:
                    for src_file in from_list:
                        source_file = str(base.joinpath(config["check"]).joinpath(src_file))
                        dest_dir = str(pub.joinpath(config["to"]))
                        shutil.copy(source_file, dest_dir)
                        print("文件 {} 拷贝至 {} 成功".format(source_file, dest_dir))
        if "function" in config.keys():
            if isinstance(config["function"], str) and config["function"][0] == '$':
                if "file" in config.keys():
                    if isinstance(config["file"], list):
                        file_path_list = [str(base.joinpath(config["check"]).joinpath(file)) for file in config["file"]]
                        if isinstance(config["new_file"], str):
                            new_file = [str(pub.joinpath(config["to"]).joinpath(config["file"]))]
                            file_path_list = file_path_list + [new_file]
                        elif isinstance(config["new_file"], list):
                            new_file_list = [str(pub.joinpath(config["to"]).joinpath(file)) for file in config["file"]]
                            file_path_list = file_path_list + new_file_list
                        globals()[config["function"][1:]](file_path_list)
                    elif isinstance(config["file"], dict):
                        file_path_dict = { param: str(base.joinpath(config["check"]).joinpath(file)) for param, file in config["file"].items()}
                        new_file_dict = { param: str(pub.joinpath(config["to"]).joinpath(file)) for param, file in config["new_file"].items()}
                        file_path_dict.update(new_file_dict)
                        globals()[config["function"][1:]](**file_path_dict)
            else:
                pass


if __name__ == "__main__":
    parse = argparse.ArgumentParser(description="根据yaml配置文件整理bsa的结果")
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
    for step, process in data.items():
        publish_files(args.indir, args.outdir, process)

