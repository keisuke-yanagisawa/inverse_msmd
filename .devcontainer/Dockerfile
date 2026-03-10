FROM kyanagis/exprorer_msmd:2025.01.25

RUN mamba install "python>=3.10" pymol-open-source -c conda-forge -y
RUN pip install rdkit libcoffee
RUN pip install mypy pytest isort ipykernel
RUN apt update
RUN apt install -y openssh-client

RUN apt install libxrender1 # for rdkit visualization

# pytest Segmentation Fault対策
# localeを適切に設定することでreadlineモジュールのクラッシュを回避
RUN apt-get update && apt-get install -y locales && \
    locale-gen en_US.UTF-8 && \
    update-locale LANG=en_US.UTF-8

# ====== settings for git / ssh =======
RUN echo "mkdir -p /root/.ssh" >> /etc/bash.bashrc
RUN echo "cp /root/.ssh_host/* /root/.ssh/" >> /etc/bash.bashrc
RUN echo "chmod 600 /root/.ssh/*" >> /etc/bash.bashrc
