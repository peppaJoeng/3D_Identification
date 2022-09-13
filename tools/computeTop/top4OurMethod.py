import pandas as pd

deli = '/'
# deli = '\\'


def parse_first(filepath, sub):
    seq = str(filepath).split(deli)
    tmp = seq[-1][:-4]
    # print(seq)
    dirc = seq[5]
    str_list = tmp.split('_')
    post = str(filepath).split(deli)[-1].split('.')[-1]
    filename = str_list[0]

    for j in range(1, len(str_list)):
        if str_list[j].isdigit():
            filename = filename + "_" + str_list[j]
            break
        filename = filename + "_" + str_list[j]
    fp = deli
    seq = str(filepath).split(deli)
    for k in range(1, len(seq) - 1):
        fp += seq[k] + deli
    fp += filename + '.' + post
    # print(fp)
    filename = fp.replace(dirc, sub)
    # print(filepath + " : " + filename)
    return filename


def parse_first4MVCNN(filepath, sub):
    seq = str(filepath).split(deli)
    #print(seq)
    tmp = seq[-1][:-4]
    # dirc = seq[5]
    dirc = seq[3]
    str_list = tmp.split('_')
    post = str(filepath).split(deli)[-1].split('.')[-1]
    filename = str_list[0]

    for j in range(1, len(str_list)):
        if str_list[j].isdigit():
            filename = filename + "_" + str_list[j]
            break
        filename = filename + "_" + str_list[j]
    fp = ""
    seq = str(filepath).split(deli)
    for k in range(0, len(seq) - 1):
        fp += seq[k] + deli
    fp += filename + '.' + post
    # print(fp)
    filename = fp.replace(dirc, sub)
    # print(filepath + " : " + filename)
    return filename


def gettopk(path, sub):
    csv = pd.read_csv(path, header=None)
    [N, M] = csv.shape
    print(str(N) + " " + str(M))
    total = 0
    total_cat = [0, 0, 0, 0, 0, 0]
    # select (topk - 1) model
    topk = [0, 2, 3, 6, 8, 11]
    ll = len(topk)

    for i in range(N):
        filename = parse_first(str(csv[0][i]), sub)
        cat = 0
        total += 1
        #print(filename)

        for k in range(1, ll):
            for j in range(topk[k - 1], topk[k]):
                if str(csv[j][i]) == filename:
                    if j == 5:
                        print(filename)
                    cat += 1
            total_cat[k - 1] += cat

    for k in range(1, ll):
        ans = (total_cat[k - 1] / total) * 100
        print(
            "top{0} = total_cat / total = {1} / {2} = {3}%".format(str(topk[k] - 1), str(total_cat[k - 1]), str(total),
                                                                   str(ans)))


if __name__ == '__main__':
    paths_corr = [
        #'result/path_crop_5_CORR_1.csv',
        #'result/path_noise_1_CORR_1.csv',
        #'result/path_quan_7_CORR_1.csv',
        #'result/path_reorder_CORR_1.csv',
        #'result/path_si_10_CORR_1.csv',
        #'result/path_sm_10_CORR_1.csv',
        'result/path_st_1_CORR_1.csv',
        #'result/path_sub_loop_CORR_1.csv',
    ]
    paths_kurt = [
        'result/path_crop_5_KURT_1.csv',
        'result/path_noise_1_KURT_1.csv',
        'result/path_quan_7_KURT_1.csv',
        'result/path_reorder_KURT_1.csv',
        'result/path_si_10_KURT_1.csv',
        'result/path_sm_10_KURT_1.csv',
        'result/path_st_1_KURT_1.csv',
        'result/path_sub_loop_KURT_1.csv',
    ]
    paths_LR = [
        'result/path_crop_5_LR_1.csv',
        'result/path_noise_1_LR_1.csv',
        'result/path_quan_7_LR_1.csv',
        'result/path_reorder_LR_1.csv',
        'result/path_si_10_LR_1.csv',
        'result/path_sm_10_LR_1.csv',
        'result/path_st_1_LR_1.csv',
        #'result/path_sub_loop_LR_1.csv',
    ]
    paths_MVCNN = [
        'MVCNN/HM25MVCNN_crop_5.csv',
        'MVCNN/HM25MVCNN_noise_1.csv',
        'MVCNN/HM25MVCNN_quan_7.csv',
        'MVCNN/HM25MVCNN_reorder.csv',
        'MVCNN/HM25MVCNN_si_10.csv',
        'MVCNN/HM25MVCNN_sm_10.csv',
        'MVCNN/HM25MVCNN_st_1.csv',
        'MVCNN/HM25MVCNN_sub_loop.csv',
    ]
    for i in range(0, len(paths_corr)):
        sub = "HM25_ori"
        # sub = "HM25MVCNN_ori"
        path = "/Users/kiki/Projects/Python/top/" + paths_corr[i]
        gettopk(path, sub)
        print(path)
