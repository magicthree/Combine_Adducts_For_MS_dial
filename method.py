import math
import os
import pandas as pd


def adduct_list(i):
    if i == 0:
        return ['[M+H]+', '[M+NH4]+', '[M+Na]+'
                # ,'[M+H-H2O]+','[M+K]+','[2M+H]+'
                ]
    if i == 1:
        return [1, 1, 1
                # ,1,1,2
                ]
    if i == 2:
        return [1.007825, 18.034374, 22.989769
                # ,-17.002739,38.963706,1.007825
                ]


def extract_file(filename, filepath, with_head):
    export_file_name = os.path.splitext(filename)[0] + "-key.csv"
    file_data = pd.read_csv(os.path.join(filepath, filename), sep='\t')
    file_data = file_data.drop(file_data.columns[file_data.apply(lambda col: 'Average' in col.values)], axis=1)
    file_data = file_data.drop(file_data.columns[file_data.apply(lambda col: 'Stdev' in col.values)], axis=1)
    if with_head:
        file_data = file_data.iloc[3:]
        file_data.columns = file_data.iloc[0]

    file_data.reset_index(drop=True, inplace=True)
    file_data = file_data.iloc[1:]
    file_data = file_data.sort_values(by=['Metabolite name', 'Adduct type', 'Average Rt(min)'])
    new_file_data = pd.concat([file_data.iloc[:, :5], file_data.iloc[:, 28], file_data.iloc[:, 32:]], axis=1)
    new_file_data = new_file_data.sort_values(by=['Metabolite name', 'Adduct type', 'Average Rt(min)'])
    new_file_data.to_csv(os.path.join(filepath, export_file_name), index=False)
    print(export_file_name + "finished")


def rm_dup(file_name, file_path, Rtol=0.1, Mtol=0.02, MtolG=0.1):
    new_file_name0 = os.path.splitext(file_name)[0] + "-filter_with_unknown.csv"
    new_file_name1 = os.path.splitext(file_name)[0] + "-filter.csv"

    file_data = pd.read_csv(os.path.join(file_path, file_name))
    file_data.reset_index(drop=True, inplace=True)

    file_data_identified = file_data.loc[
        ~file_data['Metabolite name'].str.contains('Unknown', na=False) & ~file_data['Metabolite name'].str.contains(
            'w/o', na=False)]
    file_data_identified = file_data_identified.sort_values(
        by=['Metabolite name', 'Adduct type', 'Average Rt(min)']).reset_index(drop=True)
    file_data_unknown = file_data.loc[
        file_data['Metabolite name'].str.contains('Unknown', na=False) | file_data['Metabolite name'].str.contains(
            'w/o',
            na=False)]
    file_data_unknown = file_data_unknown.sort_values(by=['Average Mz']).reset_index(drop=True)
    group = 0
    file_data_unknown.loc[:, 'group'] = -1  # 先初始化为 -1 表示还未分组

    # 手动分组，确保每个数字都与之前所有已经分组的数字进行比较
    for i in range(len(file_data_unknown)):
        if i == 0:
            # 第一个元素直接分到第 0 组
            file_data_unknown.loc[i, 'group'] = group
        else:
            # 找到与当前数字差值小于 0.1 且已经分组的行
            if any(abs(file_data_unknown.loc[i, 'Average Mz'] - file_data_unknown[file_data_unknown['group'] == group][
                'Average Mz']) < MtolG):
                file_data_unknown.loc[i, 'group'] = group
            else:
                # 否则开启新的组
                group += 1
                file_data_unknown.loc[i, 'group'] = group
    file_data_unknown = file_data_unknown.groupby('group', group_keys=False).apply(
        lambda x: x.sort_values('Average Rt(min)')).reset_index(drop=True).drop(columns=['group'])
    file_data = pd.concat([file_data_identified, file_data_unknown], ignore_index=True)

    row_count = len(file_data)
    i = 0
    row_list = []

    while i < row_count - 1:

        if file_data.loc[i, 'Metabolite name'] != file_data.loc[i + 1, 'Metabolite name']:
            row_list.append(i)
            i += 1
        elif file_data.loc[i, 'Adduct type'] != file_data.loc[i + 1, 'Adduct type']:
            row_list.append(i)
            i += 1

        elif abs(file_data.loc[i, 'Average Rt(min)'] - file_data.loc[i + 1, 'Average Rt(min)']) < Rtol and abs(
                file_data.loc[i, 'Average Mz'] - file_data.loc[i + 1, 'Average Mz']) < Mtol:
            c_line = i
            s_line = i

            while (file_data.loc[s_line, 'Metabolite name'] == file_data.loc[i + 1, 'Metabolite name'] and
                   file_data.loc[s_line, 'Adduct type'] == file_data.loc[i + 1, 'Adduct type']) and abs(
                file_data.loc[s_line, 'Average Rt(min)'] - file_data.loc[i + 1, 'Average Rt(min)']) < Rtol and abs(
                file_data.loc[i, 'Average Mz'] - file_data.loc[i + 1, 'Average Mz']) < Mtol:
                if file_data.loc[i, 'S/N average'] >= file_data.loc[i + 1, 'S/N average']:
                    c_line = i
                else:
                    c_line = i + 1
                for j in range(5, file_data.shape[1]):
                    if file_data.iloc[i, j] >= file_data.iloc[i + 1, j]:
                        file_data.iloc[i + 1, j] = file_data.iloc[i, j]
                    else:
                        file_data.iloc[i, j] = file_data.iloc[i + 1, j]
                i += 1
                if i == len(file_data) - 1:
                    break
            row_list.append(c_line)
            i += 1
        else:
            row_list.append(i)
            i += 1
    if i == row_count - 1:
        row_list.append(i)
    new_file_data = file_data.loc[row_list]
    new_file_data.to_csv(os.path.join(file_path, new_file_name0), index=False)
    new_file_data = new_file_data.loc[~(
            (file_data['Metabolite name'] == 'Unknown') | file_data['Metabolite name'].str.contains('w/o') | file_data[
        'Metabolite name'].str.contains('low score'))]
    new_file_data.to_csv(os.path.join(file_path, new_file_name1), index=False)
    print(new_file_name1 + "finished")


def alignment(filename, filepath, RTolerance=0.2):
    new_file_name = os.path.splitext(filename)[0] + "-aligned.csv"
    file_data = pd.read_csv(os.path.join(filepath, filename))
    file_data = file_data.sort_values(by=['Metabolite name', 'Adduct type', 'Average Rt(min)'])
    a_list = adduct_list(0)
    cols = ["Alignment ID", "Metabolite name", "Average Rt(min)", 'Adduct type', 'Average Mz']
    a_len = len(a_list)
    row_count = len(file_data)
    i = 0
    new_file_data = pd.DataFrame(columns=cols)
    new_rown = 0

    for i in range(0, row_count):
        if file_data.loc[i, "Metabolite name"] not in new_file_data['Metabolite name'].values:
            for adt in range(0, a_len):
                new_file_data.loc[adt + new_rown, 'Metabolite name'] = file_data.loc[i, "Metabolite name"]
                new_file_data.loc[adt + new_rown, 'Average Rt(min)'] = file_data.loc[i, "Average Rt(min)"]
                new_file_data.loc[adt + new_rown, 'Adduct type'] = a_list[adt]
                if new_file_data.loc[adt + new_rown, 'Adduct type'] == file_data.loc[i, "Adduct type"]:
                    new_file_data.loc[adt + new_rown, 'Alignment ID'] = file_data.loc[i, "Alignment ID"]
                    new_file_data.loc[adt + new_rown, 'Average Rt(min)'] = file_data.loc[i, "Average Rt(min)"]
                    new_file_data.loc[adt + new_rown, 'Average Mz'] = file_data.loc[i, "Average Mz"]
            new_rown += 1 + a_len
        else:
            namelist = new_file_data[new_file_data['Metabolite name'] == file_data.loc[i, "Metabolite name"]].index.tolist()
            already_find = False
            for j in range(0, len(namelist)):
                if abs(new_file_data.loc[namelist[j], 'Average Rt(min)'] - file_data.loc[
                    i, "Average Rt(min)"]) < RTolerance:
                    if (new_file_data.loc[namelist[j], 'Adduct type'] == file_data.loc[i, "Adduct type"]):
                        if pd.isna(new_file_data.loc[namelist[j], 'Alignment ID']):
                            new_file_data.loc[namelist[j], 'Metabolite name'] = file_data.loc[i, "Metabolite name"]
                            new_file_data.loc[namelist[j], 'Average Rt(min)'] = file_data.loc[i, "Average Rt(min)"]
                            new_file_data.loc[namelist[j], 'Average Mz'] = file_data.loc[i, "Average Mz"]
                            new_file_data.loc[namelist[j], 'Alignment ID'] = file_data.loc[i, "Alignment ID"]
                            already_find = True
                            break
            if not already_find:
                for adt in range(0, a_len):
                    new_file_data.loc[adt + new_rown, 'Metabolite name'] = file_data.loc[i, "Metabolite name"]
                    new_file_data.loc[adt + new_rown, 'Average Rt(min)'] = file_data.loc[i, "Average Rt(min)"]
                    new_file_data.loc[adt + new_rown, 'Adduct type'] = a_list[adt]
                    if new_file_data.loc[adt + new_rown, 'Adduct type'] == file_data.loc[i, "Adduct type"]:
                        new_file_data.loc[adt + new_rown, 'Alignment ID'] = file_data.loc[i, "Alignment ID"]
                        new_file_data.loc[adt + new_rown, 'Average Rt(min)'] = file_data.loc[i, "Average Rt(min)"]
                        new_file_data.loc[adt + new_rown, 'Average Mz'] = file_data.loc[i, "Average Mz"]
                new_rown += 1 + a_len

    new_file_data.to_csv(os.path.join(filepath, new_file_name), index=False)
    print(new_file_name + "finished")


def C19(filename, filepath, MTolerance=0.01, RTolerance=0.1):
    newfilename = "Standard-afterfind.csv"
    newfilename2 = "Standard-afterfind2.csv"
    filedata = pd.read_csv(os.path.join(filepath, filename))
    filedata = filedata.sort_values(by=['Metabolite name', 'Adduct type', 'Average Rt(min)'])
    Standard = [[818.665, "[M+H]+"], [524.373, "[M+H]+"], [743.709, "[M+H]+"]]
    Adductlist = adduct_list(0)
    adductM = adduct_list(1)
    Adductmz = adduct_list(2)
    Adductlen = len(Adductlist)
    cols = ["stdms", "Alignment ID", "Metabolite name", "Average Rt(min)", 'Adduct type', 'Average Mz']
    newfiledata = pd.DataFrame(columns=cols)
    newrown = 0
    for peak in Standard:
        namelist = filedata[(abs(filedata['Average Mz'] - peak[0]) < MTolerance)].index.tolist()
        MSlist = []
        for ad in range(0, len(Adductmz)):
            if peak[1] == Adductlist[ad]:
                currentad = ad

        for ad2 in range(0, len(Adductmz)):
            MSlist.append((peak[0] - Adductmz[currentad]) / adductM[ad2] * adductM[currentad] + Adductmz[ad2])

        if len(namelist) > 0:
            snlist = filedata.loc[namelist, 'S/N average']
            maxindex = snlist.idxmax()

        for adt in range(0, Adductlen):
            newfiledata.loc[adt + newrown, 'Metabolite name'] = filedata.loc[maxindex, "Metabolite name"]
            newfiledata.loc[adt + newrown, 'Average Rt(min)'] = filedata.loc[maxindex, "Average Rt(min)"]
            newfiledata.loc[adt + newrown, 'Average Mz'] = MSlist[adt]
            newfiledata.loc[adt + newrown, 'stdms'] = peak[0]
            newfiledata.loc[adt + newrown, 'Adduct type'] = Adductlist[adt]
            if newfiledata.loc[adt + newrown, 'Adduct type'] == peak[1]:
                newfiledata.loc[adt + newrown, 'Alignment ID'] = filedata.loc[maxindex, "Alignment ID"]
                newfiledata.loc[adt + newrown, 'stdms'] = peak[0]
                newfiledata.loc[adt + newrown, 'Average Rt(min)'] = filedata.loc[maxindex, "Average Rt(min)"]
                newfiledata.loc[adt + newrown, 'Average Mz'] = filedata.loc[maxindex, "Average Mz"]

        for adt in range(0, len(Adductmz)):
            if pd.isna(newfiledata.loc[adt + newrown, 'Alignment ID']):
                namelist = filedata[(abs(filedata['Average Mz'] - MSlist[adt]) < MTolerance) & (abs(
                    filedata['Average Rt(min)'] - newfiledata.loc[
                        adt + newrown, 'Average Rt(min)']) < RTolerance)].index.tolist()
                if len(namelist) > 0:
                    snlist = filedata.loc[namelist, 'S/N average']
                    maxindex = snlist.idxmax()

                    newfiledata.loc[adt + newrown, 'Metabolite name'] = filedata.loc[maxindex, "Metabolite name"]
                    newfiledata.loc[adt + newrown, 'stdms'] = peak[0]
                    newfiledata.loc[adt + newrown, 'Average Rt(min)'] = filedata.loc[maxindex, "Average Rt(min)"]
                    newfiledata.loc[adt + newrown, 'Adduct type'] = Adductlist[adt]
                    newfiledata.loc[adt + newrown, 'Alignment ID'] = filedata.loc[maxindex, "Alignment ID"]
                    newfiledata.loc[adt + newrown, 'Average Mz'] = filedata.loc[maxindex, "Average Mz"]

        newrown += 1 + Adductlen
    filedata = filedata.drop(filedata.columns[1:5], axis=1)
    newfiledata = pd.merge(newfiledata, filedata, on='Alignment ID', how='left')

    newfiledata.to_csv(os.path.join(filepath, newfilename), index=False)
    print(newfilename + "finished")

    newfiledata = pd.concat([newfiledata.iloc[:, 0], newfiledata.iloc[:, 3], newfiledata.iloc[:, 7:]], axis=1)
    new_df = pd.DataFrame(columns=newfiledata.columns)
    for i in range(0, len(newfiledata), 3):
        row_avg = newfiledata.iloc[i:i + 3, 1:2].mean()  # 前两列取平均值
        row_sum = newfiledata.iloc[i:i + 3, 2:].sum()  # 后面的列每三行取和
        new_row = pd.Series([newfiledata.iloc[i, 0]] + row_avg.tolist() + row_sum.tolist(), index=new_df.columns)
        # if row_sum[:-QCBK].sum()!=0:
        new_df = pd.concat([new_df, new_row.to_frame().T], ignore_index=True)
    new_df.to_csv(os.path.join(filepath, newfilename2), index=False)


def peakbest(name_list, t_ms, t_rt):
    if len(name_list) == 0:
        return []
    elif len(name_list) == 1:
        return [name_list.index.tolist()[0], abs(name_list.iloc[0, 0] - t_ms), abs(name_list.iloc[0, 1] - t_rt)]
    else:
        point = []
        for k in range(0, len(name_list)):
            point.append(abs(name_list.iloc[k, 0] - t_ms) + abs(name_list.iloc[k, 1] - t_rt))
        mini = point.index(min(point))
        return [name_list.index.tolist()[mini], abs(name_list.iloc[mini, 0] - t_ms),
                abs(name_list.iloc[mini, 1] - t_rt)]


def area_ratio(H, N, Na):
    lt = [H, N, Na]
    lt = [0 if x is None else x for x in lt]

    lt = [0 if math.isnan(x) else x for x in lt]

    lt = [x / max(lt) for x in lt]

    lt2 = []
    for num in lt:
        if num.is_integer():
            lt2.append(int(num))  # 如果是整数，转换为整数
        else:
            lt2.append(round(num, 3))
    return str(lt2[0]) + ":" + str(lt2[1]) + ":" + str(lt2[2])


def peak_finder(file_name, file_unknown_name, file_path, RTolerance=0.1, MTolerance=0.015):
    new_file_name = os.path.splitext(file_name)[0] + "-find.csv"
    file_data = pd.read_csv(os.path.join(file_path, file_name))

    raw_file_data = pd.read_csv(os.path.join(file_path, file_unknown_name))

    raw_file_data = raw_file_data.sort_values(by=['Metabolite name', 'Adduct type', 'Average Rt(min)'])
    raw_file_data = raw_file_data.loc[(
            (raw_file_data['Metabolite name'] == 'Unknown') | raw_file_data['Metabolite name'].str.contains('w/o') |
            raw_file_data['Metabolite name'].str.contains('low score'))]
    file_data.insert(1, "Findif", None)
    a_list = adduct_list(0)
    a_m = adduct_list(1)
    a_mz = adduct_list(2)
    a_len = len(a_list)

    cal_stat = False
    for i in range(0, len(file_data)):
        cp_n = i // a_len
        a_n = i % a_len
        std_mz = []
        if a_n == 0:
            cal_stat = False
        if pd.isna(file_data.loc[i, 'Alignment ID']):

            if not cal_stat:
                for j in range(0, a_len):
                    if not pd.isna(file_data.loc[j + cp_n * a_len, 'Alignment ID']):
                        std_mz.append(
                            (file_data.loc[j + cp_n * a_len, 'Average Mz'] - a_mz[j]) * a_m[j])
                ave_mz = sum(std_mz) / len(std_mz)

                for j in range(0, a_len):
                    if pd.isna(file_data.loc[j + cp_n * a_len, 'Alignment ID']):
                        file_data.loc[j + cp_n * a_len, 'Average Mz'] = ave_mz / a_m[j] + a_mz[j]
                cal_stat = True

            namelist = raw_file_data[(abs(raw_file_data['Average Mz'] - file_data.loc[i, "Average Mz"]) < MTolerance) & (abs(
                raw_file_data['Average Rt(min)'] - file_data.loc[i, "Average Rt(min)"]) < RTolerance)].index.tolist()

            abn = peakbest(raw_file_data.loc[namelist, ('Average Mz', 'Average Rt(min)')],
                           file_data.loc[i, "Average Mz"], file_data.loc[i, "Average Rt(min)"])
            if len(abn) != 0:
                for j in range(0, len(namelist)):
                    file_data.loc[i, "Findif"] = raw_file_data.loc[abn[0], 'Metabolite name']
                    file_data.loc[i, "Average Rt(min)"] = raw_file_data.loc[abn[0], 'Average Rt(min)']
                    file_data.loc[i, "Average Mz"] = raw_file_data.loc[abn[0], 'Average Mz']
                    file_data.loc[i, "Alignment ID"] = raw_file_data.loc[abn[0], 'Alignment ID']
                raw_file_data = raw_file_data.drop(abn[0])
    raw_file_data = pd.read_csv(os.path.join(file_path, file_unknown_name))
    raw_file_data = pd.concat([raw_file_data.iloc[:, :1], raw_file_data.iloc[:, 5:]], axis=1)
    merged_df = pd.merge(file_data, raw_file_data, on='Alignment ID', how='left')

    merged_df.to_csv(os.path.join(file_path, new_file_name), index=False)
    print(new_file_name + "finished")


def blank_filter(file_name, file_path, min, ratio, qc, bk_list):
    file_data = pd.read_csv(os.path.join(file_path, file_name))
    file_data = file_data[[col for col in file_data.columns if col not in qc] + qc]
    file_data = file_data[[col for col in file_data.columns if col not in bk_list] + bk_list]
    bk_len = len(bk_list)

    new_file_name = os.path.splitext(file_name)[0] + "-bk.csv"
    for index, row in file_data.iterrows():
        if not pd.isna(file_data.loc[index, 'Alignment ID']):
            max_value = row[bk_list].max() * ratio
            if max_value < min:
                max_value = min
            file_data.iloc[index, 7:-bk_len] = row[7:-bk_len].where(row[7:-bk_len] > max_value, 0)
            file_data.to_csv(os.path.join(file_path, new_file_name), index=False)
    print(new_file_name + "finished")


def merge_file(file_name, file_path, qc, bk_list):
    a_len = len(adduct_list(0))
    qcbk = len(qc) + len(bk_list)
    new_file_name = os.path.splitext(file_name)[0] + "-mged.csv"
    file_data = pd.read_csv(os.path.join(file_path, file_name))
    a_data = file_data.iloc[:, 0:6]
    print(file_data)

    file_data = pd.concat(
        [file_data.iloc[:, 2], file_data.iloc[:, 4], file_data.iloc[:, 5], file_data.iloc[:, 3], file_data.iloc[:, 7:]],
        axis=1)
    print(file_data)
    new_df = pd.DataFrame(columns=file_data.columns)
    for i in range(0, len(file_data), a_len):
        row_avg = file_data.iloc[i:i + a_len, 3:4].mean()  # 前两列取平均值
        row_sum = file_data.iloc[i:i + a_len, 4:].sum()  # 后面的列每三行取和
        adduct = ""
        mass = ""
        for j in range(i, i + a_len):
            if pd.notna(a_data.iloc[j, 0]):
                adduct = adduct + a_data.iloc[j, 4] + ";"
                mass = mass + str(a_data.iloc[j, 5]) + ";"
        adduct = adduct[:-1]
        mass = mass[:-1]
        new_row = pd.Series([file_data.iloc[i, 0]] + [adduct] + [mass] + row_avg.tolist() + row_sum.tolist(),
                            index=new_df.columns)

        if row_sum[:-qcbk].sum() != 0 or qcbk == 0:
            new_df = pd.concat([new_df, new_row.to_frame().T], ignore_index=True)
    new_df = new_df.rename(columns={'Average Mz': 'Mz'})

    new_df.to_csv(os.path.join(file_path, new_file_name), index=False)
    print(new_file_name + "finished")
