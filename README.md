# HMM-Predict-Startcodon
使用HMM和bootstrap预测酵母起始密码子，注释基因组

1、相关数据准备
Step1中实验物种Candida albicans（SC5341）的基因组序列文件(.fna)和注释文件(.gff)。
从EPD数据库中下载TATA-box HMM数据，选取表头为’TATA-box HMM trained from 900 unrelated general promoter sequences’表格，复制粘贴保存为TATA-boxHMM.txt文件。

2、DNA 元件的计算鉴别
编写TATA鉴别的python程序step5.py：
以TATA-boxHMM.txt文件数据创建TATA-box的HMM模型，默认wordsize为12bp，步长为1bp，p值cutoff0.05，bootstrap中打乱次数为100次，并进行次数为10，cutoff为0.5的预处理，对基因组中未测出序列(N)，非明确碱基(W S等)，均忽略不计入在内，保留score，-log(score)，p值，染色体序号，起始位点，正负链，与预测的序列，输出为out.gff文件。


将相关文件（.py .fna hmm模型）打包压缩为zip文件上传至小型机
scp E:\step5.zip s31@42.244.7.51:
mv step5.zip genomic/
unzip step5.zip
cd step5.zip
nohup python3 step5.py & 
 
再编写step5.py并使用1000行的基因组序列进行测试时，发现其结果的数量随着实验次数，存在无规律的波动性，其原因是bootstrap算法本身存在随机性，12bp4种碱基的组成非常多，只打乱1000次，再加上碱基组成的偏差，只能算作部分抽样，故存在随机性。
修改cutoff与阈值的参数，重复上述实验，对参数及结果记录如下表（仅记录实验成功的数据）
表格 4 tata预测参数及结果数据表
cutoff	打乱次数	结果数量
0.05	100	282710
0.05	1000	297524
0.01	200	80922
0.005	100	77429
0.005	500	52268
0.005	1000	44617
#由基因组注释文档可知，该物种已知基因共6245个
