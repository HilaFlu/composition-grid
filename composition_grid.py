import numpy as np
import pandas as pd
import itertools
import argparse
import os
import sys
import time

my_parser = argparse.ArgumentParser(prog='composition_grid', description='Creates updated composition grid files. '
                                                                         'for 0 concentration type in file: "-1 -1 1,'
                                                                         'to keep concentration type in file: 0 0 1"',
                                    epilog="Good luck with your experiments! :)")
my_parser.add_argument('Data_Path', metavar='data_path', type=str, help='the path to the lipodimics data')
my_parser.add_argument('Setup_Path', metavar='setup_path', type=str, help='the path to setup data')
my_parser.add_argument('New_Grid_Path', metavar='new_grid_path', type=str, help="the path to the new grid's file")
args = my_parser.parse_args()
new_grid_path = args.New_Grid_Path
old_grid_path = args.Data_Path
wanted_lipids_file_path = args.Setup_Path

if not os.path.isdir(new_grid_path):
    print('The output path specified does not exist: %s' % new_grid_path)
    sys.exit()

if not os.path.isfile(old_grid_path):
    print("Couldn't find input grid data in file: %s" % old_grid_path)
    sys.exit()

if not os.path.isfile(wanted_lipids_file_path):
    print("Couldn't find config file: %s" % wanted_lipids_file_path)
    sys.exit()

"""new_grid_path = "path/files"
old_grid_path = "/path/lipidomics_with_chol_scaled.csv"
wanted_lipids_file_path = "/path/config_file.csv""""


def set_concentration_factor(x, membrane_part, new_con, old_con):
    # finds the new concentration factor for one wanted lipid
    total_sum = np.sum(x[membrane_part])
    print("new con:",new_con, "old con:", old_con, "total:", total_sum)
    factor = (total_sum - ((total_sum / 100) * new_con)) / (total_sum - old_con)
    print(factor, membrane_part)
    return factor


def calc_old_con(x, changed_lipids, membrane_part):
    # calculates the old concentration for all wanted lipids in layer and returns their sum
    con, index_num = 0, 0
    if membrane_part == "inside":
        index_num = 1
    elif membrane_part == "outside":
        index_num = 2
    zero_index_list = []
    for item_lipid in changed_lipids:
        for item_name in item_lipid[0]:
            value = x.loc[x['ID'] == item_name, [membrane_part]]
            con = con + np.sum(value[membrane_part])
            if item_lipid[index_num] != [0.0, 0.0, 1.0] or item_lipid[index_num][0] >= 0:
                zero_index_list.append(changed_lipids.index(item_lipid))
    return con, zero_index_list


def calc_new_con_grid(x, membrane_part, wanted_lipids_list):
    index_num = 0
    if membrane_part == "inside":
        index_num = 1
    elif membrane_part == "outside":
        index_num = 2
    new = []
    for lipid in wanted_lipids_list:
        if lipid[index_num][0] < 0:
            old_con_keep = 0
            total_sum = np.sum(x[membrane_part])
            for item_name in lipid[0]:
                old_con_keep = old_con_keep + (x.loc[x['ID'] == item_name, [membrane_part]])
            if index_num == 1:
                lipid[index_num] = [(old_con_keep.inside.item()/total_sum)*100]
            else:
                lipid[index_num] = [(old_con_keep.outside.item()/total_sum)*100]
        else:
            lipid[index_num] = np.arange(lipid[index_num][0], lipid[index_num][1] + lipid[index_num][2],
                                         lipid[index_num][2])
        new.append(lipid[index_num])
    array = []
    first, *rest = new
    for tup in itertools.product(first, *rest):
        array.append(tup)
    return array


def process_config_file(wanted_lipids_file):
    global old_grid
    wanted_lipids_list = wanted_lipids_file.values.tolist()
    working_lipids = []
    selected_lipids_list = []
    print(wanted_lipids_list)
    for item in wanted_lipids_list:
        if item[0] != item[0]:
            continue
        item[0] = item[0].split()
        if item[4] == item[4]:
            item[4] = item[4].split()
        if item[3] == "manual":
            if item[2] != item[2] or item[1] != item[1]:
                raise Exception("Please specify wanted concentration variables for both leaflets :)")
            item[2] = [float(s) for s in item[2].split()]
            item[1] = [float(s) for s in item[1].split()]
        elif item[3] == "fixed":
            if (item[1] == item[1]) and (item[2] == item[2]):
                raise Exception("Please specify wanted concentration variables only for the chosen leaflet :)")
            out_grid = 0.0
            in_grid = 0.0
            for lipid in item[0]:
                out_grid = out_grid + ((old_grid.loc[old_grid["ID"] == str(lipid),
                                                     ["outside"]].values[0][0]) / np.sum(old_grid["outside"])) * 100
                print(out_grid, type(out_grid))
                in_grid = in_grid + ((old_grid.loc[old_grid["ID"] == str(lipid),
                                                   ["inside"]].values[0][0]) / np.sum(old_grid["inside"])) * 100
                print(in_grid, out_grid, "hj")
            if item[1] == item[1]:
                non_listed_layer = 2
                listed_layer = 1
                listed_grid = in_grid
                non_listed_grid = out_grid
            else:
                non_listed_layer = 1
                listed_layer = 2
                listed_grid = out_grid
                non_listed_grid = in_grid
            item[listed_layer] = [float(s) for s in item[listed_layer].split()]
            if (out_grid > 0) and (out_grid is not np.nan) and (in_grid is not np.nan) \
                    and (in_grid > 0):
                ratio = non_listed_grid / listed_grid
                item[non_listed_layer] = [round(float(s) * ratio, 8) for s in item[listed_layer]]
            elif (non_listed_grid == 0) or (non_listed_grid is np.nan):
                item[non_listed_layer] = [0.0, 0.0, 1.0]
            elif (listed_grid is np.nan) or (listed_grid == 0):
                item[non_listed_layer] = [float(s) for s in item[listed_layer]]
                item[listed_layer] = [0.0, 0.0, 1.0]
        print(item[4])
        if type(item[4]) == float:
            if item[1] == [0.0, 0.0, 1.0] or item[2] == [0.0, 0.0, 1.0] or item[1][0] < 0 or item[2][0] < 0:
                working_lipids.insert(0, [item[0], item[1], item[2]])
            else:
                working_lipids.append([item[0], item[1], item[2]])
        else:
            if item[1] == [0.0, 0.0, 1.0] or item[2] == [0.0, 0.0, 1.0] or item[1][0] < 0 or item[2][0] < 0:
                selected_lipids_list.insert(0, [item[0], item[1], item[2], item[4]])
            else:
                selected_lipids_list.append([item[0], item[1], item[2], item[4]])
    return working_lipids, selected_lipids_list


def set_concentration_value(x, changed_lipid, new_con, membrane_part, total_sum, factor):
    # calculates and returns the new concentration for the lipid
    if factor == 0:
        factor = 1
    """index_layer = 0
    if membrane_part == "inside":
        index_layer = 1
    elif membrane_part == "outside":
        index_layer = 2"""
    ids = []
    for lipid in changed_lipid:
        if (len(lipid[0]) > 1) and (str(x['ID']) in lipid[0]):
            ids.append(lipid[0])
            index_num = ids.index(lipid[0])
            if new_con[index_num] >= 0:
                new_concentration = calc_fam_factor(lipid[0], new_con[index_num], membrane_part, total_sum,
                                                    str(x['ID']))
            else:
                new_concentration = x[membrane_part] * factor
            return new_concentration
        else:
            ids.append(lipid[0][0])
    if str(x['ID']) in ids:
        index_num = ids.index(x['ID'])
        if new_con[index_num] >= 0:
            new_concentration = (total_sum / 100) * new_con[index_num]
        else:
            new_concentration = x[membrane_part] * factor
    else:
        new_concentration = x[membrane_part] * factor
    return new_concentration


def calc_fam_factor(lipids, concentration, membrane_part, total_sum, chosen_lipid):
    # calculates and returns the family concentration factor
    fam_sum = 0
    for lipid in lipids:
        fam_sum = fam_sum + old_grid.loc[old_grid["ID"] == str(lipid), [membrane_part]].values[0][0]
    chosen_con = old_grid.loc[old_grid["ID"] == chosen_lipid, [membrane_part]].values[0][0]
    fam_factor = ((total_sum * (concentration / 100)) / fam_sum) * chosen_con
    return fam_factor


def set_selected_mode_value(x, membrane_part, new_con, lipids_selection):
    total_sum = np.sum(new_grid[membrane_part])
    # create ID list for [0] and [3]
    ids_0 = []
    ids_3 = []
    for lipid in lipids_selection:
        if str(x['ID']) in lipid[0]:
            ids_0.append(lipid[0])
            index_num = ids_0.index(lipid[0])
            check_if_selected(membrane_part, new_con[index_num], lipid)
            if new_con[index_num] >= 0:
                original_part_sum_after = (total_sum / 100) * new_con[index_num]
                original_part_sum_before = 0
                for original_lipid in lipid[0]:
                    original_part_sum_before = original_part_sum_before + \
                                               (new_grid.loc[new_grid["ID"] == str(original_lipid),
                                                             [membrane_part]].values[0][0])
                original_con = (new_grid.loc[new_grid["ID"] == str(x["ID"]), [membrane_part]].values[0][0])
                new_selected_concentration = (original_part_sum_after / original_part_sum_before) * original_con
            else:
                new_selected_concentration = x[membrane_part]
            return new_selected_concentration
        elif str(x['ID']) in lipid[3]:
            ids_3.append(lipid[3])
            index_num = ids_3.index(lipid[3])
            check_if_selected(membrane_part, new_con[index_num], lipid)
            if new_con[index_num] >= 0:
                original_part_sum_after = (total_sum / 100) * new_con[index_num]
                selected_part_sum_before = 0
                for selected_lipid in lipid[3]:
                    selected_part_sum_before = selected_part_sum_before + \
                                               (new_grid.loc[new_grid["ID"] == str(selected_lipid),
                                                             [membrane_part]].values[0][0])
                original_part_sum_before = 0
                for original_lipid in lipid[0]:
                    original_part_sum_before = original_part_sum_before + \
                                               (new_grid.loc[new_grid["ID"] == str(original_lipid),
                                                             [membrane_part]].values[0][0])
                selected_con = (new_grid.loc[new_grid["ID"] == str(x["ID"]), [membrane_part]].values[0][0])
                # print(selected_con, original_part_sum_before, original_part_sum_after, selected_part_sum_before, "mmmm")
                new_selected_concentration = (((original_part_sum_before + selected_part_sum_before) -
                                               original_part_sum_after) / selected_part_sum_before) * selected_con
            else:
                new_selected_concentration = x[membrane_part]
            return new_selected_concentration
        else:
            new_selected_concentration = x[membrane_part]
            ids_0.append(lipid[0][0])
            ids_3.append(lipid[3][0])
    return new_selected_concentration


def check_if_selected(membrane_part, new_con, selected):
    total_sum = np.sum(old_grid[membrane_part])
    primer_sum = 0
    for primer_lipid in selected[0]:
        primer_sum = primer_sum + (old_grid.loc[old_grid["ID"] == str(primer_lipid), [membrane_part]].values[0][0])
    new_con_for = (total_sum / 100) * new_con
    select_sum = 0
    for lipid in selected[3]:
        select_sum = select_sum + (old_grid.loc[old_grid["ID"] == str(lipid), [membrane_part]].values[0][0])
    if ((select_sum < (new_con_for - primer_sum)) and (new_con_for > primer_sum)) or \
            ((primer_sum < (new_con_for - select_sum)) and (new_con_for < primer_sum)):
        print(select_sum, primer_sum, new_con_for)
        raise Exception("Selected lipids concentration is too small :(")


def match_in_out(inner_array, outer_array, inner_index, outer_index):
    matched_con = []
    max_match = len(inner_array) / len(outer_array)
    print(inner_array, "kfg", outer_array, "match")
    in_small = False
    if max_match < 1:
        max_match = 1 / max_match
        in_small = True
    if in_small:
        small_array, small_index, small_leaflet = inner_array, inner_index, 0
        big_array, big_index, big_leaflet = outer_array, outer_index, 1
    else:
        small_array, small_index, small_leaflet = outer_array, outer_index, 1
        big_array, big_index, big_leaflet = inner_array, inner_index, 0
    for little_con in small_array:
        insert = [[], []]
        insert[small_leaflet] = little_con
        matched_con.append(insert)
    for con_slot in matched_con:
        match_count = max_match
        regulator = 0
        trash_can = [big_array[0]]
        print(small_index, "small")
        while match_count > 0:
            big_con = big_array[regulator]
            bool_val = []
            for bigger_index in small_index:
                if big_con[bigger_index] == trash_can[-1][bigger_index]:
                    bool_val.append(True)
                else:
                    bool_val.append(False)
            if all(bool_val):
                match_count -= 1
                con_slot[big_leaflet].append(big_con)
                trash_can.append(big_array.pop(regulator))
                print(trash_can[-1], big_con, regulator, match_count, bool_val)
                regulator -= 1
            regulator += 1
            print(big_array, regulator, small_array, trash_can[-1])
            if len(big_array) == 0:
                break
    return matched_con, big_leaflet


# def positive_sum(concentrations):
#     pos_sum = 0
#     for item in concentrations:
#         print(item, "checkup")
#         if item > 0:
#             pos_sum += item
#     return pos_sum


# set updated composition grid file
# old_grid = pd.read_csv("/Users/flomin/Desktop/composition_grid/setup_files/lipidomics_with_chol.csv")
old_grid = pd.read_csv(old_grid_path)
old_grid = old_grid.fillna(0)
# wanted_lipid_file = pd.read_csv("/Users/flomin/Desktop/composition_grid/setup_files/config_file.csv")
wanted_lipid_file = pd.read_csv(wanted_lipids_file_path)
# membrane_parts = input("Choose membrane layer you want to change. \nFor both layers type 'total_concentration'."
#                       "\nFor outer layer type 'outside'. \nFor inner layer type 'inside':")
# wanted_lipids = [[lipid_name, inside_tuple, outside_tuple], [...], [...]]
wanted_lipids, selected_lipids = process_config_file(wanted_lipid_file)
print(wanted_lipids, "sep", selected_lipids)
folder_name = str(time.strftime("%Y-%m-%d--%H-%M-%S"))
folder_path = os.path.join(new_grid_path, folder_name)
os.mkdir(folder_path)
# wanted_lipids = [["Cer4412", (3, 5, 1), (3, 5, 1)], ["Cer4202", (2, 2.5, 0.5), (0, 0, 1)],
#                 ["Cer4222", (1, 2, 1), (2, 4, 2)]]
follow = 0
count = 0
if len(selected_lipids) > 0:
    blank1, index_out_select = calc_old_con(old_grid, selected_lipids, "outside")
    blank2, index_in_select = calc_old_con(old_grid, selected_lipids, "inside")
    array_out_selected = calc_new_con_grid(old_grid, "outside", selected_lipids)
    array_in_selected = calc_new_con_grid(old_grid, "inside", selected_lipids)
    matched_array_selected, big_leaf = match_in_out(array_in_selected, array_out_selected,
                                                    index_in_select, index_out_select)
if len(wanted_lipids) > 0:
    old_con_out, index_out = calc_old_con(old_grid, wanted_lipids, "outside")
    old_con_in, index_in = calc_old_con(old_grid, wanted_lipids, "inside")
    array_out = calc_new_con_grid(old_grid, "outside", wanted_lipids)
    array_in = calc_new_con_grid(old_grid, "inside", wanted_lipids)
    # print(array_out, "-", array_in)
    # print(array_out_selected, "---", array_in_selected)
    # print(selected_lipids)
    # print(wanted_lipids[0], index_in, index_out)
    # print(matched_array_selected, big_leaf, "matched")
    if len(array_in) >= len(array_out):
        big = array_in
        small = array_out
        chosen = "in"
        zero_index_small = index_out
        zero_index_big = index_in
    else:
        big = array_out
        small = array_in
        chosen = "out"
        zero_index_big = index_out
        zero_index_small = index_in
    mod_factor = len(big) / len(small)
    print(len(big), "  ", len(small))
    for small_con in small:
        trash = [big[0]]
        current_mod_factor = mod_factor
        reg = 0
        while current_mod_factor > 0:
            bigger_con = big[reg]
            print("reg: " + str(reg) + " last: " + str(bigger_con) + "mod: " + str(current_mod_factor), "small: ",
                  str(small_con))
            val_bool = []
            for index_big in zero_index_small:
                if bigger_con[index_big] == trash[-1][index_big]:
                    val_bool.append(True)
                else:
                    val_bool.append(False)
            if all(val_bool):
                count += 1
                current_mod_factor -= 1
                if chosen == "in":
                    inner_con = bigger_con
                    outer_con = small_con
                else:
                    inner_con = small_con
                    outer_con = bigger_con
                # new_grid = pd.read_csv("/Users/flomin/Desktop/composition_grid/setup_files/lipidomics_with_chol.csv")
                new_grid = pd.read_csv(old_grid_path)
                if len(wanted_lipids) > 0:
                    new_grid["outside"] = old_grid.apply(set_concentration_value, changed_lipid=wanted_lipids,
                                                         new_con=outer_con,
                                                         membrane_part="outside",
                                                         total_sum=np.sum(old_grid["outside"]),
                                                         factor=set_concentration_factor(
                                                             old_grid, "outside", sum(outer_con), old_con_out),
                                                         axis="columns")
                    new_grid["inside"] = old_grid.apply(set_concentration_value, changed_lipid=wanted_lipids,
                                                        new_con=inner_con,
                                                        membrane_part="inside",
                                                        total_sum=np.sum(old_grid["inside"]),
                                                        factor=set_concentration_factor(
                                                            old_grid, "inside", sum(inner_con), old_con_in),
                                                        axis="columns")
                    new_grid["total_concentration"] = \
                        new_grid["inside"].fillna(value=0) + new_grid["outside"].fillna(value=0)
                    if len(selected_lipids) == 0:
                        # new_grid.to_csv(folder_path + "/new_composition_grid" + wanted_lipids[0][0][0] + "_in" +
                        #                str(sum(inner_con)) + "_out" + str(sum(outer_con)) + "_" + str(count) + ".csv")
                        new_grid.to_csv(os.path.join(folder_path, str(follow) + ".csv"))
                        follow += 1
                        trash.append(big.pop(reg))
                        print(trash[-1])
                        reg -= 1
                        break
                    new_selected_grid = new_grid.copy()
                    for tot_concentration in matched_array_selected:
                        for big_concentration in tot_concentration[big_leaf]:
                            follow += 1
                            if big_leaf == 1:
                                new_selected_grid["inside"] = new_grid.apply(set_selected_mode_value,
                                                                             membrane_part="inside",
                                                                             new_con=tot_concentration[0],
                                                                             lipids_selection=selected_lipids,
                                                                             axis="columns")
                                new_selected_grid["outside"] = new_grid.apply(set_selected_mode_value,
                                                                              membrane_part="outside",
                                                                              new_con=big_concentration,
                                                                              lipids_selection=selected_lipids,
                                                                              axis="columns")
                            else:
                                new_selected_grid["inside"] = new_grid.apply(set_selected_mode_value,
                                                                             membrane_part="inside",
                                                                             new_con=big_concentration,
                                                                             lipids_selection=selected_lipids,
                                                                             axis="columns")
                                new_selected_grid["outside"] = new_grid.apply(set_selected_mode_value,
                                                                              membrane_part="outside",
                                                                              new_con=tot_concentration[1],
                                                                              lipids_selection=selected_lipids,
                                                                              axis="columns")
                            new_selected_grid["total_concentration"] = (new_selected_grid["inside"].fillna(value=0)
                                                                        + new_selected_grid["outside"].fillna(value=0))
                            # new_selected_grid.to_csv(folder_path + "/new_composition_grid" +
                            #                          wanted_lipids[0][0][0] + "_in" + str(sum(inner_con)) + "_out" +
                            #                          str(sum(outer_con)) + "_" + str(follow) + ".csv")
                            new_selected_grid.to_csv(os.path.join(folder_path, str(follow) + ".csv"))
                trash.append(big.pop(reg))
                print(trash[-1])
                reg -= 1
            reg += 1
            if len(big) == 0:
                break
    print(follow)
elif len(selected_lipids) > 0:
    new_grid = pd.read_csv(old_grid_path)
    new_selected_grid = new_grid.copy()
    for tot_concentration in matched_array_selected:
        for big_concentration in tot_concentration[big_leaf]:
            follow += 1
            if big_leaf == 1:
                new_selected_grid["inside"] = new_grid.apply(set_selected_mode_value,
                                                             membrane_part="inside",
                                                             new_con=tot_concentration[0],
                                                             lipids_selection=selected_lipids,
                                                             axis="columns")
                new_selected_grid["outside"] = new_grid.apply(set_selected_mode_value,
                                                              membrane_part="outside",
                                                              new_con=big_concentration,
                                                              lipids_selection=selected_lipids,
                                                              axis="columns")
            else:
                new_selected_grid["inside"] = new_grid.apply(set_selected_mode_value,
                                                             membrane_part="inside",
                                                             new_con=big_concentration,
                                                             lipids_selection=selected_lipids,
                                                             axis="columns")
                new_selected_grid["outside"] = new_grid.apply(set_selected_mode_value,
                                                              membrane_part="outside",
                                                              new_con=tot_concentration[1],
                                                              lipids_selection=selected_lipids,
                                                              axis="columns")
            new_selected_grid["total_concentration"] = (new_selected_grid["inside"].fillna(value=0)
                                                        + new_selected_grid["outside"].fillna(value=0))
            # new_selected_grid.to_csv(folder_path + "/new_composition_grid" +
            #                          wanted_lipids[0][0][0] + "_in" + str(sum(inner_con)) + "_out" +
            #                          str(sum(outer_con)) + "_" + str(follow) + ".csv")
            new_selected_grid.to_csv(os.path.join(folder_path, str(follow) + ".csv"))
    print(follow)
