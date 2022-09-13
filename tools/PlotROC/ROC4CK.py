import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import auc


def Positive(path, thred):
    csv = pd.read_csv(path, header=None)
    [N, M] = csv.shape
    # print(str(N) + " " + str(M))
    if N != M:
        print("ERROR : Matrix is not N * N")
        return 0, 0
    TP = 0
    FN = 0
    for i in range(N):
        if float(csv[i][i]) >= 1000:
            continue
        # >= thred <==> true
        # < thred <==> false
        # real positive and predicted positive (TP)
        if float(csv[i][i]) >= thred:
            TP = TP + 1
        # real positive and predicted negative (FN)
        else:
            FN = FN + 1
    return TP, FN


def Negative(path, thred):
    csv = pd.read_csv(path, header=None)
    [N, M] = csv.shape
    # print(str(N) + " " + str(M))
    if N != M:
        print("ERROR : Matrix is not N * N")
        return 0, 0
    FP = 0
    TN = 0
    for i in range(N):
        for j in range(N):
            if float(csv[i][j]) >= 1000:
                continue
            elif i == j:
                continue
            # real negative and predicted positive (FP)
            elif float(csv[i][j]) >= thred:
                FP = FP + 1
            # real negative and predicted negative (TN)
            else:
                TN = TN + 1
    return FP, TN


def ROC(path):
    TPR = []
    FNR = []
    FPR = []
    TNR = []
    Thred = []
    for thred in np.arange(0.9, 1, 0.0001):
        Thred.append(thred)
        # print(thred)
        [TP, FN] = Positive(path, thred)
        TPR.append(TP / (TP + FN))
        FNR.append(FN / (TP + FN))
        [FP, TN] = Negative(path, thred)
        FPR.append(FP / (FP + TN))
        TNR.append(TN / (FP + TN))
    return TPR, FPR, Thred


def ROC2(path):
    TPs = []
    FNs = []
    FPs = []
    TNs = []
    Threds = []
    for thred in np.arange(0.9, 1, 0.0001):
        Threds.append(thred)
        # print(thred)
        [TP, FN] = Positive(path, thred)
        TPs.append(TP)
        FNs.append(FN) 
        [FP, TN] = Negative(path, thred)
        FPs.append(FP)
        TNs.append(TN)

    # print(TPR)
    # print(FPR)
    return TPs, FNs, FPs, TNs, Threds


def plotROC(Thred, TPR, FPR, title='Characteristic'):
    plt.scatter(x=Thred, y=TPR, label='TPR', color='r')
    plt.scatter(x=Thred, y=FPR, label='FPR', color='g')
    plt.legend(loc='lower right')

    plt.title(title)
    # plt.plot([(0, 0), (1, 1)], 'r--')
    plt.xlim([0.8999, 1])
    plt.ylim([0, 1.01])
    plt.ylabel('Rate')
    plt.xlabel('Thresholds')

    mm = 0
    ind = 0
    for i in range(len(FPR)):
        if (TPR[i] - FPR[i] > mm):
            mm = TPR[i] - FPR[i]
            ind = i
    print(Thred[ind])
    print(TPR[ind])
    print(FPR[ind])
    plt.plot([Thred[ind], Thred[ind]], [TPR[ind], FPR[ind]], 'b-')
    plt.show()


def plotfigure(Threds, TPRs, FPRs, titles):
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    l = len(TPRs)
    r = 2
    c = l / 2
    for i in range(l):
        ax = fig.add_subplot(r, c, i + 1)
        ax.scatter(x=Threds[i], y=TPRs[i], label='TPR', color='r')
        ax.scatter(x=Threds[i], y=FPRs[i], label='FPR', color='g')
        ax.legend(loc='lower right')

        # ax.title(titles[i])
        # plt.plot([(0, 0), (1, 1)], 'r--')
        plt.xlim([0.8999, 1])
        plt.ylim([0, 1.01])
        plt.ylabel('Rate')
        plt.xlabel('Thresholds')
        # plt.show()
    fig.show()


def titleName(path):
    # print(path.split('/')[1])
    # print(path.split('/')[1].split('.')[0])
    filename = path.split('/')[1].split('.')[0].split('_')
    # print(filename)
    ans = ''
    for i in range(2, len(filename)):
        ans += filename[i]
    return ans


def singleRate(distance_LR):
    TPRs = []
    FPRs = []
    Threds = []
    titles = []
    for i in range(0, len(distance_LR)):
        path = "/Users/kiki/Projects/Python/top/" + distance_LR[i]

        [TPR, FPR, Thred] = ROC(path)
        TPRs.append(TPR)
        FPRs.append(FPR)
        Threds.append(Thred)
        titles.append(titleName(distance_LR[i]))
        # [TP, FP] = Positive(path, 0.01)
        # print("the true positive is : " + str(TP))
        # print("the false positive is : " + str(FP))
    return TPRs, FPRs, Threds, titles


def allRate(distance_LR):
    ATPs = []
    AFNs = []
    AFPs = []
    ATNs = []
    Threds = []
    for i in range(0, len(distance_LR)):
        path = "/Users/kiki/Projects/Python/top/" + distance_LR[i]
        [TPs, FNs, FPs, TNs, Thred] = ROC2(path)
        ATPs.append(TPs)
        AFNs.append(FNs)
        AFPs.append(FPs)
        ATNs.append(TNs)
        Threds.append(Thred)
    TP = np.sum(ATPs, axis=0)
    FN = np.sum(AFNs, axis=0)
    FP = np.sum(AFPs, axis=0)
    TN = np.sum(ATNs, axis=0)

    SP = np.sum([TP, FN], axis=0)
    SN = np.sum([FP, TN], axis=0)

    TPR = [a / b for a, b in zip(TP, SP)]
    FPR = [a / b for a, b in zip(FP, SN)]

    return TPR, FPR, Threds[0]


if __name__ == '__main__':

    distance_CORR = [
        'CORR/distance_HM25_crop5_CORR_1.csv',
        'CORR/distance_HM25_noise_1_CORR_1.csv',
        'CORR/distance_HM25_quan_7_CORR_1.csv',
        'CORR/distance_HM25_reorder_CORR_1.csv',
        'CORR/distance_HM25_si_10_CORR_1.csv',
        'CORR/distance_HM25_sm_10_CORR_1.csv',
        'CORR/distance_HM25_st1_CORR_1.csv',
    ]
    # distance_KURT = [
    #     'KURT/distance_HM25_crop5_KURT_1.csv',
    #     'KURT/distance_HM25_noise_1_KURT_1.csv',
    #     'KURT/distance_HM25_quan_7_KURT_1.csv',
    #     'KURT/distance_HM25_reorder_KURT_1.csv',
    #     'KURT/distance_HM25_si_10_KURT_1.csv',
    #     'KURT/distance_HM25_sm_10_KURT_1.csv',
    #     'KURT/distance_HM25_st1_KURT_1.csv',
    # ]

    [TPR, FPR, Thred] = allRate(distance_CORR)
    # [Threds, TPRs, FPRsm, titles] = singleRate(distance_LR)

    # plotfigure(Thred=Threds, TPRs=TPRs, FPRs=FPRs, titles=titles)
    plotROC(Thred, TPR, FPR)
