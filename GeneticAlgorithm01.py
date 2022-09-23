# coding:utf-8

"""
editor: hap
time: 2022年9月22日
history:happp-hap/improve
"""

"""

问题：

公式：y = a*x**3 + b*x**2 + c*x + d

在闭区间[-1,3]上，用遗传算法求最大值。其中
a=0.27, b=-0.4, c=-0.86, d=2 
  

问题输入：
xx条xx位二进制“遗传基因”

每条xx位二进制“遗传基因”，转换成十进制后，代表x坐标轴上的x坐标。


问题输出：
01条xx位二进制“遗传基因”

代表最优解的x值，输出转换后的十进制数字，并输出极大值点。

"""


import time
import random
import logging
import math
import os
import matplotlib.pyplot as plt
import numpy as np
from logging import handlers

class Logger(object):
    """
    配置日志模块
    """
    level_relations = {
        'debug':logging.DEBUG,
        'info':logging.INFO,
        'warning':logging.WARNING,
        'error':logging.ERROR,
        'crit':logging.CRITICAL
    }#日志级别关系映射

    def __init__(self, filename, level='info', when='D', backCount=3, \
        fmt='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'):

        self.logger = logging.getLogger(filename)

        format_str = logging.Formatter(fmt)#设置日志格式

        self.logger.setLevel(self.level_relations.get(level))#设置日志级别

        sh = logging.StreamHandler()#往屏幕上输出
        sh.setFormatter(format_str) #设置屏幕上显示的格式

        #往文件里写入#指定间隔时间自动生成文件的处理器
        th = handlers.TimedRotatingFileHandler(filename=filename, when=when, backupCount=backCount, encoding='utf-8')
        #实例化TimedRotatingFileHandler
        #interval是时间间隔，backupCount是备份文件的个数，如果超过这个个数，就会自动删除，when是间隔的时间单位，单位有以下几种：
        # S 秒
        # M 分
        # H 小时、
        # D 天、
        # W 每星期（interval==0时代表星期一）
        # midnight 每天凌晨
        #设置文件里写入的格式
        th.setFormatter(format_str)

        #把对象加到logger里s
        self.logger.addHandler(sh) 
        self.logger.addHandler(th)

if __name__ == '__main__':
    log = Logger('impove.log',level='warning')
    #log.logger.debug('debug')
    #log.logger.info('info')
    #log.logger.warning('警告')
    #log.logger.error('报错')
    #log.logger.critical('严重')
    #Logger('error.log', level='error').logger.error('error')


class Improve():
    """
    遗传基因进化类
    """

    # f QUESTION
    A  = 0.27
    B  = -0.4
    C  = -0.86
    D  = 2
    # 求解问题的x轴左右区间
    X0 = -1
    X1 = 3

    # 基因种类,letter=1、digit=2、binary=3
    GENE_TYPE             = 3
    # 优化精度，两代之间小于优化精度，结束进化
    OPTIMIZATION_ACCURACY = 1e-6
    # x轴取值粒度 Granularity 1e-8时，优化精度1e-6，在导数低于100的情况下可能达到
    X_GRANULARITY         = 1e-8
    # 环境承载量
    # ENV_LOAD              = 20
    # 种群数量, 设置表示固定种群数量, 本参数目前需要双数
    GROUP_NUM             = 80
    # 进化代数上限（达到上限不再进化）
    IMPROVE_MAX_EPOCH     = 500
    #基因长度
    GEN_LEN               = int(math.log((X1-X0)/X_GRANULARITY, 2)+1);
    # 杂交概率 0.4~0.9
    CROSSOVER_PROBABILITY = 0.5
    # 变异概率 0.0001~0.1
    CROSSOVER_PROBABILITY = 0.01
    # 进化速度
    # IMPROVE_SPEED       = 200
    # 上帝的恩惠
    # GOD_GRACE           = 3
    # 环境承载压力补偿 0.1~1
    # ENV_PRESSURE_OFFSET = 0.5
    # ENV_CREATION
    # ENV_CREATION        = 0.999
    # 保存文件的位置
    SAVE_FILE_PATH        = "saves.txt"


    def createGen(self, length, type_=3):
        """
        产生一个基因序列
        如：1000011000110100100010
        """

        def getALetter():
            """随机产生一位letter的闭包"""
            letter = "abcdefghijklmnopqrstuvwxyz"
            k = int(random.random()*len(letter))
            return letter[k]

        def getADigit():
            """随机产生一位Digit的闭包"""
            digit  = "0123456789"
            k = int(random.random()*len(digit))
            return digit[k]

        def getABinary():
            """随机产生一位二进制的闭包"""
            binary  = "01"
            k = int(random.random()*len(binary))
            return binary[k]

        # 根据type_选择基因类型
        type_switch = {1:getALetter, 2:getADigit, 3:getABinary}
        getAElement = type_switch[type_]

        # 基因是随机产生的，随机产生每一位，共产生length位
        gen = ""
        for i in range(length):
            gen += getAElement()
        log.logger.info("生成一个随机gen，为："+gen)
        return gen

    def createGeneToFile(self, number):
        gens = []
        f = open(self.SAVE_FILE_PATH, "w+")

        # 生成 GROUP_NUM 个长度为 GEN_LEN 的基因
        for i in range(number):
            
            while True:
                gen = self.createGen(self.GEN_LEN, self.GENE_TYPE)
                x = self.f_map(self.binary2decimal(gen))

                # 在x上，二进制上下界，不见得完全匹配十进制上下界，
                # 这里排除掉界外的值，重新生成
                if self.X0 < x < self.X1:
                    break

            gens.append(gen)

        self.writeGenToFile(gens)
        log.logger.info("生成"+str(len(gens))+"个个体")
        return gens

    def readGenFromFile(self):
        """
        从文件中读取现存的基因列表,若文件不存在，则创建
        """
        gens = []

        if os.path.exists(self.SAVE_FILE_PATH):
            f = open(self.SAVE_FILE_PATH, "r")
            
            for each in f.readlines():
                if each != "" or each != "\n" or each != "\r\n" or each != "\r":
                    if each.strip() != "":
                        gens.append(each.strip())
            f.close()
            if len(gens) == 0:
                gens = self.createGeneToFile(self.GROUP_NUM)

            else:
                log.logger.info('从'+str(self.SAVE_FILE_PATH)+"中读取"+str(len(gens))+"个个体")

        else:
            gens = self.createGeneToFile(self.GROUP_NUM)

        
        return gens

    def writeGenToFile(self, gens):
        """将基因序列保存到文件"""
        # f = open(self.SAVE_FILE_PATH, "w")
        # for each in gens:
        #     f.write(each+"\n")
        # f.close()
        write_gens = []
        for each in gens:
            write_gens.append(each+"\n")
        f = open(self.SAVE_FILE_PATH, "w")
        f.writelines(write_gens)
        f.close()

    def binary2decimal(self, binary):
        """
        输入str类型二进制
        输出数字十进制
        """
        decimal = int(binary,2)*1e-6
        return decimal

    def decimal2binary(self,decimal):
        """
        输出str类型的二进制
        """
        binary = str(bin(int(int(decimal)/self.X_GRANULARITY)))[2:]
        return binary

    def f(self,binary_input):
        decimal = self.binary2decimal(binary_input)
        x = self.f_map(decimal)
        y = self.A*x**3 + self.B*x**2 + self.C*x + self.D
        return y

    def f_map(self,decimal):
        return decimal-1

    def normalized(self, y_list):
        sum_ = sum(y_list)
        log.logger.debug("标准化前："+str(y_list))
        if sum_ == 0 :
            y_list_normalized = [ 0 for y in y_list]
        else:
            y_list_normalized = [ y/sum_ for y in y_list]
        log.logger.debug("标准化后："+str(y_list_normalized))
        return y_list_normalized

    def roulette(self, gens):

        y_list = [ self.f(x) for x in gens]
        min_ = min(y_list)

        # 轮盘赌算法需要标准化位正数
        if min_<0 :
            # 标准化只需要正数，上移
            # +up_offset,否则造成最小值淘汰，防止最小值无被选择的概率
            up_offset = abs((max(y_list) - min(y_list))/len(y_list))
            y_list = [ y-min_+up_offset for y in y_list] 

        log.logger.debug("标准化")
        y_list = self.normalized(y_list)

        y_list_roulette = []
        sum_ = 0
        for each in y_list:
            sum_ += each
            y_list_roulette.append(sum_)

        rand_ = random.random()
        
        i = 0
        for each in y_list_roulette:

            if each > rand_:
                break;
            i+=1

        return gens[i]

    def select(self, gens):

        assert(len(gens)>1)

        new_gens = []
        for i in range(self.GROUP_NUM):
            new_gens.append(self.roulette(gens))

        return new_gens

    def crossover(self, gens):
        """
        交换
        1和11，2和12，3和13，以此类推，10和20
        1. 随机位置
        2. 多个位置 todo if needed
        """
        def exchange(gens, i, j):

            k = int(random.random()*self.GEN_LEN)
            gen1 = gens[i][:k] + gens[j][k:]
            gen2 = gens[i][k:] + gens[j][:k]

            # 在x上，二进制上下界，不见得完全匹配十进制上下界，
            # 这里排除掉界外的值，使用原本未杂交的基因
            gen1_x = self.f_map(self.binary2decimal(gen1))
            gen2_x = self.f_map(self.binary2decimal(gen2))
            if  self.X0 < gen1_x < self.X1: 
                gens[i] = gen1
            if  self.X0 < gen2_x < self.X1:
                gens[j] = gen2

        for i in range(int(self.GROUP_NUM/2)):
            if random.random() > self.CROSSOVER_PROBABILITY :
                exchange(gens, i, i+10)
        return gens

    def mutation(self, gens):
        """
        变异
        """
        j=0
        for each in gens:
            for i in range(len(each)):
                if (random.random() < self.CROSSOVER_PROBABILITY):
                    if each[i]=="1":
                        each = each[:i]+"0"+each[i+1:]
                    else:
                        each = each[:i]+"1"+each[i+1:]
                    log.logger.info("第"+str(j)+"条基因，第"+str(i)+"位已发生突变！")
            j+=1
        return gens

    def get_best(self, gens):
        y_list = [self.f(x) for x in gens]
        max_ = max(y_list)

        k = 0
        for each in y_list:
            if(each==max_):
                break
            k+=1
        return gens[k]

    def viewAllGens(self,gens):
        for each in gens:
         log.logger.info(str(each)+ "  x="+str(self.f_map(self.binary2decimal(each))) + "  y="+str(self.f(each)))

    def plot_f_func(self):
        x = np.arange(-1, 3, 0.1) # x轴采样
        y = 0.27*x**3-0.4*x**2-0.86*x+2       # 计算对应的 y
        plt.title('y=0.27*x^3-0.4*x^2-0.86*x+2')    # 图像名称
        plt.xlabel('x')      # x轴标签
        plt.ylabel('y')      # y轴标签
        plt.plot(x, y)     # 绘制图像
        plt.show()        # 显示图像

    def improve(self,):
        """
        种群迭代器
        """
        epoch = 1;
        best_last = float('-inf')
        stop = 10


        # read_from_file 读取
        # 从文件中读取现有的 gens, 没有则随机生成
        gens = self.readGenFromFile()
        log.logger.info("种群大小："+str(self.GROUP_NUM))
        log.logger.info("基因长度："+str(self.GEN_LEN))
        self.viewAllGens(gens)
        while True :
            
            # select 选择
            log.logger.info("选择")
            gens = self.select(gens)
            

            self.viewAllGens(gens)

            # crossover 交换
            log.logger.info("交换")
            gens = self.crossover(gens)
            
            self.viewAllGens(gens)

            # mutation 变异
            log.logger.info("变异")
            self.mutation(gens)
            
            self.viewAllGens(gens)

            # optimization 调优

            gen = self.get_best(gens)
            best_slution = self.f(gen)
            optimization = best_slution - best_last
            best_last = best_slution
           
            log.logger.warning("当前是第"+ str(epoch) + "代:\n最优基因为"+ str(gen) +"，\
                x = "+  str(self.f_map(self.binary2decimal(gen))) +" best slution "+str(best_slution)+" optimization:"+str(optimization))

            # 本次代没有产生新的突破，等待5代
            if abs(optimization) == 0:
                stop-=1
            else:
                stop = 10

            if(abs(optimization) < self.OPTIMIZATION_ACCURACY):
                if stop < 0:
                    log.logger.warning("小于预设精度"+ str(self.OPTIMIZATION_ACCURACY) + "，不再迭代，最优化完成")
                    return 

            if(epoch > self.IMPROVE_MAX_EPOCH):
                log.logger.warning("代数超过最大值"+ str(self.IMPROVE_MAX_EPOCH) + "，未达到最优化，停止迭代")
                return 

            epoch += 1

        self.writeGenToFile(gens)


if __name__ == '__main__':
    Improve().improve()
    Improve().plot_f_func()