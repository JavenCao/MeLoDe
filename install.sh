#!/bin/bash

SCRIPT_DIR=$(realpath $(dirname "$0"))
echo "Script absolute directory: $SCRIPT_DIR"

# 检查是否安装了 Go
if command -v go &> /dev/null; then
    echo "Go is already installed. Version: $(go version)"
    #exit 0
else
    echo "Go is not installed. Proceeding with installation."
    # 设置 Go 的版本
    GO_VERSION=1.20.4

    # 设置安装路径
    INSTALL_DIR="$SCRIPT_DIR/mygo"

    # 设置下载目录
    DOWNLOAD_DIR="$SCRIPT_DIR/Downloads"
    mkdir -p $DOWNLOAD_DIR

    # 下载 Go 安装包
    wget -q -O $DOWNLOAD_DIR/go$GO_VERSION.tar.gz https://go.dev/dl/go$GO_VERSION.linux-amd64.tar.gz

    # 解压 Go 安装包
    mkdir -p $INSTALL_DIR
    tar -C $INSTALL_DIR -xzf $DOWNLOAD_DIR/go$GO_VERSION.tar.gz --strip-components=1


    # 创建工作目录
    mkdir -p $SCRIPT_DIR/go_workspace

    # 配置环境变量
    #echo "export PATH=\$PATH:$INSTALL_DIR/bin" >> ~/.bashrc
    #echo "export GOPATH=\$HOME/go_workspace" >> ~/.bashrc
    #echo "export GOROOT=$INSTALL_DIR" >> ~/.bashrc
    export PATH=$PATH:$INSTALL_DIR/bin
    export GOPATH=$SCRIPT_DIR/go_workspace
    export GOROOT=$INSTALL_DIR

    # 使配置生效
    # source ~/.bashrc

    # 删除下载的安装包
    #rm $DOWNLOAD_DIR/go$GO_VERSION.tar.gz

fi


export CGO_CFLAGS="-I$SCRIPT_DIR/libtensorflow/include"
echo $CGO_CFLAGS

#decompress libtensorflow.so.2.13.0.tar.bz2
tar -xvjf $SCRIPT_DIR/libtensorflow/lib/libtensorflow.so.2.13.0.tar.bz2 -C $SCRIPT_DIR/libtensorflow/lib/

export LIBRARY_PATH=$SCRIPT_DIR/libtensorflow/lib:$LIBRARY_PATH
echo $LIBRARY_PATH

#go mod tidy
bash -c "cd ./src/read_predict && go mod tidy"
bash -c "cd ./src/site_predict && go mod tidy"
bash -c "cd ./src/tag_decode && go mod tidy"


go build -C ./src/read_predict -o ../../bin/MeLoDe_Read cmd/cmd.go

go build -C ./src/site_predict -o ../../bin/MeLoDe_Combo cmd/combo.go

go build -C ./src/site_predict -o ../../bin/MeLoDe_HiFi cmd/hifi.go

go build -C ./src/tag_decode -o ../../bin/MeLoDe_TagExtract cmd/cmd.go