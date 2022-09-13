import pandas as pd

deli = '\\'


def parse_first(filepath):
    seq = str(filepath).split(deli)
    tmp = seq[-1][:-4]
    str_list = tmp.split('_')
    post = str(filepath).split(deli)[-1].split('.')[-1]
    filename = str_list[0]

    for j in range(1, len(str_list)):
        if str_list[j].isdigit():
            filename = filename + "_" + str_list[j]
            break
        filename = filename + "_" + str_list[j]
    fp = 'D:' + deli
    seq = str(filepath).split(deli)
    for k in range(1, len(seq) - 1):
        fp += seq[k] + deli
    fp += filename + '.' + post
    # print(fp)
    filename = fp
    # print(filepath + " : " + filename)
    return filename


def gettopk(path):
    csv = pd.read_csv(path, header=None)
    [N, M] = csv.shape
    print(str(N) + " " + str(M))
    total = 0
    total_cat = [0, 0, 0, 0, 0, 0]
    # select (topk - 1) model
    topk = [0, 2, 3, 6, 8, 11]
    ll = len(topk)

    for i in range(N):
        filename = parse_first(str(csv[0][i]))
        cat = 0
        total += 1
        #print(filename)

        for k in range(1, ll):
            for j in range(topk[k - 1], topk[k]):
                if str(csv[j][i]) == filename:
                    cat += 1
            total_cat[k - 1] += cat

    for k in range(1, ll):
        ans = (total_cat[k - 1] / total) * 100
        print(
            "top{0} = total_cat / total = {1} / {2} = {3}%".format(str(topk[k] - 1), str(total_cat[k - 1]), str(total),
                                                                   str(ans)))


if __name__ == '__main__':
    paths_meshnet = [
        'Meshnet/crop_5.csv',
        'Meshnet/noise_1.csv',
        'Meshnet/quan_7.csv',
        'Meshnet/reorder.csv',
        'Meshnet/si_10.csv',
        'Meshnet/sm_10.csv',
        'Meshnet/st_1.csv',
        'Meshnet/sub_loop.csv',
    ]
    paths_USC = [
        'USC/crop_5.csv',
        'USC/noise_1.csv',
        'USC/quan_7.csv',
        'USC/reorder.csv',
        'USC/si_10.csv',
        'USC/sm_10.csv',
        'USC/st_1.csv',
        'USC/sub_loop.csv',
    ]
    for i in range(0, len(paths_USC)):
        path = "/Users/kiki/Projects/Python/top/" + paths_USC[i]
        gettopk(path)
        print(path)
