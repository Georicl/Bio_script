import subprocess
import sys
import argparse
import os
#############################
#这是用R画图的一个脚本，可以画火山图与热图与PCA，并且可以选择安装orgDB
#############################
###region 脚本建立
#获取根目录
project_root = os.path.dirname(os.path.abspath(__file__))
#获取脚本解释器
draw_plot = argparse.ArgumentParser(
    description = "这个脚本用于转录组处理结果的处理，需要搭配已有的RNAseq.sh使用后，将得出的结果的3.merge和4.DE_result文件进行处理，最后会导出为PDF文件，文件位置为当前目录。"
)
###endregion
###region 参数设定
#版本信息
draw_plot.add_argument('--version',action='version',version='version: 1.0')
#添加参数
draw_plot.add_argument('-m','--merge',type=str,help='请输入merge_result文件夹路径')
draw_plot.add_argument('-o','--org',help='请输入orgDB文件路径，需要orgDB源文件，默认为无')
draw_plot.add_argument('-S','--SAMPLE',help='显示是否将表达量改名为组名，选择改名请输入文件路径,默认为无')
draw_plot.add_argument('-d','--de',help='请输入DE_result文件')
draw_plot.add_argument('-s','--sample',help='请输入sample.txt文件，该文件要求只有两列，第一列为组名，第二列为样本名，不要标题')
#参数解析
read_args = draw_plot.parse_args()
# 如果没有任何参数，输出--help
if len(sys.argv) == 1:
    draw_plot.print_help()
    sys.exit(1)  # 退出程序，返回状态码1表示错误
#R脚本路径读取
path_stv = os.path.join(project_root,'draw_plot_R','save_the_result.R')
path_org = os.path.join(project_root,'draw_plot_R','make_orgDB.R')
path_draw = os.path.join(project_root, 'draw_plot_R', 'draw_plot_R.R')
path_pca = os.path.join(project_root, 'draw_plot_R','PCA.R')
### endregion
### region进行脚本处理
#如果需要对exp表达量改名
if read_args.SAMPLE:
    name_exp_sample = "1"
else:
    name_exp_sample = "0"
#如果有orgDB，则触发处理orgDB的脚本
if read_args.org:
    subprocess.run(['Rscript',path_org,read_args.org])
    choose_org = "1"
else:
    choose_org = "0"
#运行save_the_result.R
subprocess.run(['Rscript',path_stv,read_args.merge,read_args.de,name_exp_sample,read_args.SAMPLE,path_org,choose_org])
#输出画图文件
subprocess.run(['Rscript', path_draw])
#输出PCA聚类
subprocess.run(['Rscript',path_pca,read_args.sample])
### endregion
#捕获当前目录
current_directory = os.getcwd()
#删除文件
file_delete = os.path.join(current_directory, "RNAseq.rdata")
#检测文件是否存在
if os.path.exists(file_delete):
    os.remove(file_delete)

#产生日志
#将输出重定向到文件 rscript_output.log
with open('draw_plot_R.log', 'w') as f:
    sys.stdout = f  #捕获标准输出
    sys.stderr = f  #捕获错误输出
