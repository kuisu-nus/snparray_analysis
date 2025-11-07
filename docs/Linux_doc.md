# Linux 记录

## Docker 
### Install

```bash
#!/bin/bash
echo "安装 Docker NVIDIA 支持..."

# 检查系统
distribution=$(. /etc/os-release;echo $ID$VERSION_ID)
echo "检测到系统: $distribution"

# 安装 NVIDIA Container Toolkit
echo "1. 安装 NVIDIA Container Toolkit..."
curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg
curl -s -L https://nvidia.github.io/libnvidia-container/$distribution/libnvidia-container.list | \
    sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \
    sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list

sudo apt update
sudo apt install -y nvidia-container-toolkit

# 配置 Docker
echo "2. 配置 Docker..."
sudo nvidia-ctk runtime configure --runtime=docker

# 重启服务
echo "3. 重启 Docker 服务..."
sudo systemctl restart docker

# 验证
echo "4. 验证安装..."
sudo docker run --rm --gpus all nvidia/cuda:12.0-base nvidia-smi

echo "安装完成！"
```

- 精简版本测试
```bash
sudo docker run --rm --gpus all nvidia/cuda:12.0.0-runtime-ubuntu20.04 nvidia-smi
```

- 数据挂在到docker
```bash
docker run -it --rm \
    -v /path/to/your/project/code:/app/code \
    -v /path/to/your/data:/app/data \
    your_image_name:tag \
    /bin/bash

# --rm: 容器停止后自动移除
# --it: 交互式运行，并分配一个伪终端
# /bin/bash: 启动容器后直接进入Bash Shell

# 假设您在 /Users/user/my_plink_project 目录下
docker run -it --rm \
    -v $(pwd):/workspace/code \
    your_plink_image:latest
# 此时，宿主机上的 /Users/user/my_plink_project 映射到容器内的 /app

docker run -it --rm \
    -v /home/sukui/03.projects/02.gene/Genos/:/workspace/code/Genos/ \
    zhanxiaoai/mega:v1 \
    /bin/bash

# use GPU
docker run -it --rm --gpus all \
    -v /home/sukui/03.projects/02.gene/Genos/:/workspace/code/Genos/ \
    zhanxiaoai/mega:v1 \
    /bin/bash
```

### 操作容器
- `docker ps -a`: 查看所有的容器
- `docker stop`: 正常停止一个正在运行的容器
- `docker kill`: 强制终止一个容器
- `docker pause`: 暂停容器中的所有进程
- `docker exec -it`: 在正在运行的容器中执行命令，常用于获取shell
- `docker log -f`: 查看容器的实时标准输出日志, 类似于tail -f
- `docker  inspect`: 查看容器的详细配置
- `docker stats`: 实时监控容器的CPU，内存等
- `docker rm -f`: 强制移除正在运行的容器


## 文本文件处理

- 计数
```bash
# 统计列的分布
awk '{count[NF]++} END{for(i in count) print i, count[i]}' filename | sort -n
# 统计指定列的分布
awk '{count[$5]++} END{for(i in count) print i, count[i]}' PGT_TLS_ALL1.refactor_step8.ped | sort -n
```

- 文件处理
```bash
# 只保留VCF文件的表头和数据行
awk '/^#CHROM/ || !/^#/' input.vcf > output.vcf

# 提取指定的列
awk '{print $3}' input.vcf > snp_ids.txt

awk -F',' 'NR==1 {print "SNP,OR,P"} NR>1 && $9 < 0.05 {print $2 "," $7 "," $9}' file.csv > significant.csv
# `NR==1`: NR表示行号，这里匹配第一行
# `{print "SNP,OR,P"}`: 输出新的表头"SNP, OR, P"
# `NR>1`: 匹配第2行及以后的数据行
# `$9<0.05`: 条件判断，第9列<0.05
 # `{print $2 "," $7 "," $9}`：输出第2列(SNP)、第7列(OR)、第9列(P)，用逗号连接
```

## PLINK操作
- 从vcf中提取指定的SNPs
```bash
# 基于snp_list.txt
plink --vcf input.vcf --extract snplist.txt --recode vcf --out output
# 基于position.txt
plink --vcf input.vcf --extract range position.txt --recode vcf --out output


```
