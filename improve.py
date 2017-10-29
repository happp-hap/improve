# coding:utf-8
import time
import random
import logging

def logCfg():
    logging.basicConfig(level=logging.DEBUG,
                format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                datefmt='%a, %d %b %Y %H:%M:%S',
                filename='improve.log',
                filemode='w')
    # logging.debug('This is debug message')
    # logging.info('This is info message')
    # logging.warning('This is warning message')

class Improve():

    LETTER_TYPE         = 1
    DIGIT_TYPE          = 2
    ENV_LOAD            = 50
    GEN_LEN             = 10
    IMPROVE_SPEED       = 200
    GOD_GRACE           = 3
    ENV_PRESSURE_OFFSET = 0.5
    ENV_CREATION        = 0.999
    SAVE_FILE_PATH      = "saves.txt"

    def readGenFromFile(self):
        """
        从文件中读取现存的基因列表
        """
        f = open(self.SAVE_FILE_PATH, "r")
        gens = []
        for each in f.readlines():
            if each != "" or each != "\n" :
                gens.append(each[:len(each)-1])
        f.close()
        return gens

    def createGen(self, type_, length):
        """产生一个基因"""

        def getALetter():
            """随机产生letter的闭包"""
            letter = "abcdefghijklmnopqrstuvwxyz"
            k = int(random.random()*len(letter))
            return letter[k]

        def getADigit():
            """随机产生Digit的闭包"""
            digit  = "0123456789"
            k = int(random.random()*len(digit))
            return digit[k]

        # 根据type_选择基因类型
        type_switch = {1:getALetter, 2:getADigit}
        getAElement = type_switch[type_]

        # 基因是随机产生的，随机产生每一位，共产生length位
        gen = ""
        for i in range(length):
            gen += getAElement()
        return gen

    def getParents(self, gens):
        """
        从种群中获得双亲
        1. 可能存在自交
        2. 禁止自交
        """
        assert(len(gens)>1)

        k1 = int(random.random()*len(gens))
        k2 = int(random.random()*len(gens))
        if (k1 != k2):
            genFath = gens[k1]
            genMon = gens[k2]
        else:
            k2 = (k1 + 1) % len(gens)
            genFath = gens[k1]
            genMon = gens[k2]
        return genFath, genMon

    def genRecombination(self, genFath, genMon):
        """
        重组
        """
        assert((str == type(genFath)) and (self.GEN_LEN == len(genFath))
             and (self.GEN_LEN == len(genMon)) and (str == type(genMon)))

        # k = int(random.random()*self.GEN_LEN)
        k = int(self.GEN_LEN/2)
        gen = genFath[:k] + genMon[k:]
        return gen

        
    def envFeatures(self, it, nowGensNum):
        """
        做特征的判断和存活概率的生成
        此环境的特征选择方式为：
        1. 相同。即判断相邻的基因是否相同。相同则存活概率大。
        2. 相同。基因最后一个元素和第一个元素做相同比较。
        """

        # 断言基因的类型
        assert(str == type(it))

        # 环境特征的判断
        same_num = 0
        for i in range(len(it)):
            if it[i%len(it)] == it[(i+1)%len(it)] :
                same_num += 1

        # 存活概率的生成
        # 存活概率 = (基因中相同的个数/基因长度) * 
        # (1-(现存个体数/环境承载)*环境承载压力补偿 ) 
        prob = ( same_num/len(it)*self.GOD_GRACE ) * ( 1-(nowGensNum/self.ENV_LOAD)*self.ENV_PRESSURE_OFFSET ) 

        return prob

    def randomEventRun(self, prob, times):
        """
        随机时间运行器
        input: 概率，运行次数
        """
        return bool( (random.random()*times) < prob )

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

    def improve(self,):
        """
        种群迭代器
        """
        while True :
            ## BEGIN： 生成gens-种群的全部成员           ##
            # 从文件中读取现有的 gens
            gens = self.readGenFromFile()
            nowGensNum = len(gens)
            logging.info('从'+self.SAVE_FILE_PATH+"中读取到"+str(nowGensNum)+"个个体")

            # 产生一个gen
            # 字母类型，GEN_LEN位
            if (random.random() < self.ENV_CREATION or len(gens)<2):
                gen = self.createGen(self.LETTER_TYPE, self.GEN_LEN)
                logging.info("生成一个随机gen，为："+gen)
            else:
                genFath, genMon = self.getParents(gens)
                gen = self.genRecombination(genFath, genMon)
                logging.info("生成一个重组gen，为："+gen)         

            # 将gen加到gens中，现在gens是整个种群的gens
            # 种群每次进化产生一个个体
            gens.append(gen)
            nowGensNum += 1
            logging.info("将新生成的gen加到种群中，现在种群大小："+str(nowGensNum))
            ## END  ： 生成gens-种群的全部成员           ##
            
            # 条件判断
            # 执行随机环境裁决，传入概率和裁决次数，返回1为生
            next_gens = []
            for each in gens :
                prob = self.envFeatures(each, nowGensNum)
                # 活下来的组成下一个种群
                if (self.randomEventRun(prob, 1)):
                    next_gens.append(each)

            # 种群进化需要时间
            logging.info("种群进化需要时间: "+str(1/self.IMPROVE_SPEED))
            time.sleep(1/self.IMPROVE_SPEED)

            # 存活的新的种群
            self.writeGenToFile(next_gens)
            logging.info("存活的新的种群"+str(len(next_gens))+"个")

if __name__ == '__main__':
    logCfg()
    Improve().improve()
