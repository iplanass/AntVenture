# Tested on python3.11-3.13
# Investigate distribution of food & foraging dynamics depending on foraging method (T vs SB vs both), distance bw nest, terrain, food source (viscosity) & colony size (or prop of foragers)
# Author: Isaac Planas-Sitja
# Last modification: 2024/03/07
# Name:AntVenture: Sims

# All parameters have been estimated using real data from Fujioka H, Marchand M, and LeBoeuf AC. "Diacamma ants adjust liquid foraging strategies in response to biophysical constraints." Proceedings of the Royal Society B 290.2000 (2023): 20230549.

import numpy as np
import random
import pandas as pd
import time
import warnings
import matplotlib.pyplot as plt
# ~ import PySimpleGUI as sg
import os
from datetime import datetime
from scipy import stats
# ~ import _tkinter
# ~ from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

warnings.filterwarnings('ignore')

#Define model
def diacamma_model(N, Nf, D, time_sim, visco='NA', sugar='NA', behavior=1, terrain=0, method_sb="simple", n_sims=1):
    # TIME AND VOLUME
    t_sb = 10  # time to grab social bucket 
    tt_t = 75  # time to grab trophallaxis
    v_sb = 10  # volume for social bucket ; v_t depends on tt_t and viscosity
    fed_vol = 10  # quantity to transfer in order to feed 1 ant
    # speed_t = 1  # walking speed trophallaxis
    # PROBABILITIES
    p_out = 0.1  # probability to go out and start exploring
    p_nest = 0.1  # probability go back to nest
    p_source = 0.2  # probability find food source
    p_feed_t = 0.1  # p_feed_sb not used, as it is directly computed in the model
    warning_message = "No errors found"
    if behavior not in range(0, 3):
        warning_message = "Behaviour is not defined."  # -- add n 2 = first tropha, and then grab, both beh at same time
        return None, warning_message
    Nn = N - Nf  # Number of nurses
    foragers = np.arange(0, Nf)
    # ---------------------------#
    #	Terrain modification	#
    # ---------------------------#
    p_drop_sb = 0.001 * (1 + terrain)
    speed_sb = 1 / (1 + terrain)
    dist_t = int(D)  # time needed to travel bw nest - source = distance/speed -- no need, speed_t = 1
    dist_sb = int(D / speed_sb)  # time needed to travel bw nest - source = distance/speed
    dist_both = int(D / speed_sb)  # use same speed as trophallaxis or SB, we can choose
    # -------------------------------#
    #	Volume & time modification	#
    # -------------------------------#
    if visco != 'NA':
        t_t = int(tt_t - 1.4 * visco * 1000)
        v_t = round((1 / (31.377 + 7043.306 * visco)) * t_t * 10, 1)  # multiply by 10 because the model works in 10muL
        #
        t_both = t_sb + t_t
        v_both = v_sb + v_t
        if v_t < 0 or t_t < 0:
            warning_message = "Viscosity (mPA) is too high, please choose a lower value."
            return None, warning_message
    elif sugar != 'NA':
        visco = (0.3074 * np.exp(7.29 * sugar)) / 1000  # divided to get PA
        t_t = int(tt_t - 1.4 * visco * 1000)
        v_t = round((1 / (31.377 + 7043.306 * visco)) * t_t * 10, 1)  # multiply by 10 because the model works in 10muL
        t_both = t_sb + t_t
        v_both = v_sb + v_t
        if v_t < 0 or t_t < 0:
            warning_message = "Sugar concentration is too high, please choose a lower value [0-1]."
            return None, warning_message
    else:
        warning_message = "Please introduce a value for either Viscosity (mPA) or Sugar (sugar concentration in range 0-1)."
        return None, warning_message
    # -----------------------#
    #	START: CREATE TABLE	#
    # -----------------------#
    feature_list = ["task", "behav", "where", "p_temp", "timing", "max_time", "fed", "naif", "qliquid"]
    dat = pd.DataFrame(0, index=np.arange(N), columns=feature_list)
    dat.loc[foragers, "task"] = 0  # forager
    dat.loc[Nf:, "task"] = 1  # nurse
    dat.loc[:, "behav"] = behavior  # 1=tropha; 0=SB
    results = np.zeros([int(time_sim / 10 + 1), 7])  # needs to be the same as vector below
    results[0] = [0, N, 0, 0, 0, 0, 0]
    vec_time = np.arange(0, time_sim + 1, 10)  # we do +1 to include 100
    t0 = time.time()
    fin_res = []
    for l in np.arange(0, n_sims):
        k = 0  # timing
        results[:] = 0
        results[:, 6] = int(l)
        dat.iloc[:, 2:] = 0
        for t in np.arange(1, time_sim + 1):
            # SELECT NAIF INDIVIDUALS DEPENDING ON WHERE
            Nf_in = np.where(dat.loc[foragers, ["where", "naif"]].sum(axis=1) == 0)[0]  # who is in nest and naif
            Nf_out = np.where((dat.loc[foragers, "where"] == 1) & (dat.loc[foragers, "naif"] == 0))[
                0]  # who is outside and naif
            Nf_source = np.where(dat.loc[foragers, "where"] == 2)[0]  # who is at the source
            Nf_ret = np.where(dat.loc[foragers, "where"] == 3)[0]  # who is going back to nest
            Nf_tropha = np.where(
                (dat.loc[foragers, "where"] == 0) & (dat.loc[foragers, "naif"] == 1) & (dat.loc[foragers, "behav"] == 1) & (
                        dat.loc[foragers, "qliquid"] > 0))[0]  # find foragers with liquid to share troph
            Nf_sb = np.where(
                (dat.loc[foragers, "where"] == 0) & (dat.loc[foragers, "naif"] == 1) & (dat.loc[foragers, "behav"] == 0) & (
                        dat.loc[foragers, "qliquid"] > 0))[0]  # find foragers with liquid to share SB
            Nf_trsb = np.where(
                (dat.loc[foragers, "where"] == 0) & (dat.loc[foragers, "naif"] == 1) & (dat.loc[foragers, "behav"] == 2) & (
                        dat.loc[foragers, "qliquid"] > 0))[0]  # find foragers with liquid to share tropha + SB
            # -------------------------#
            # 1- ANTS INSIDE NEST phase#
            # -------------------------#
            # 1.1 see if naif go out
            if len(Nf_in) > 0:
                dat.loc[Nf_in, "p_temp"] = np.random.random(len(Nf_in)) - p_out  # random prob for all naives inside nest
                dat.loc[Nf_in, "where"] = np.where((dat.loc[Nf_in, "p_temp"] < 0), 1,
                                                   dat.loc[Nf_in, "where"])  # those that are negative, go outside
            # 1.2. informed empty ants go out
            temp = np.where(
                (dat.loc[foragers, "where"] == 0) & (dat.loc[foragers, "naif"] == 1) & (dat.loc[foragers, "qliquid"] == 0))[
                0]  # who is in nest, Informed (naif=1) & empty --> goes out
            for i in temp:
                if dat.loc[i, "behav"] == 1:  # if tropha
                    dat.loc[i, ["where", "timing"]] = [1,
                                                       -1]  # we put -1 because in the phase "outside", they will get +1 timing, so they will get at zero and the end of timestep
                elif dat.loc[i, "behav"] == 0:  # SB
                    dat.loc[i, ["where", "timing"]] = [1, -1]  # max time already in dist_sb/t mode
                else:  # both
                    dat.loc[i, ["where", "timing"]] = [1, -1]  # max time already in dist_sb/t mode]
            dat.loc[:, "p_temp"] = 0  # reboot proba to 0
            # -------------------------------------------#
            # 2- TROPHALLAXIS & SOCIAL BUCKET INSIDE NEST#
            # -------------------------------------------#
            ######################
            # 2.1 - Trophallaxis  #
            ######################
            if len(Nf_tropha) > 0:
                Nempty = np.where((dat.loc[:, "fed"] < 1) & (dat.loc[:, "where"] == 0))[0]
                if len(Nempty) > 0:  # if empty ants in nest
                    dat.loc[Nf_tropha, "p_temp"] = np.random.random(len(Nf_tropha)) - p_feed_t * len(
                        Nempty) / N  # 1/feeding_time * probability to find empty ant; could be N/(area*something else). Makes up for time passing food --------------- check max_capacity
                    if len(np.where((dat.loc[Nf_tropha, "p_temp"] < 0) & (dat.loc[Nf_tropha, "qliquid"] >= fed_vol))[
                               0]) > 0:
                        temp = random.sample(sorted(Nempty), len(
                            (np.where((dat.loc[Nf_tropha, "p_temp"] < 0) & (dat.loc[Nf_tropha, "qliquid"] >= fed_vol)))[
                                0]))  # we choose randomly which individual is fed with fed_vol
                        dat.loc[temp, "fed"] = 1  # those are fed
                        temp = np.where((dat.loc[Nf_tropha, "p_temp"] < 0) & (dat.loc[Nf_tropha, "qliquid"] >= fed_vol))
                        dat.loc[Nf_tropha[temp], "p_temp"] = 0
                        dat.loc[Nf_tropha[temp], "qliquid"] = dat.loc[Nf_tropha[temp], "qliquid"] - fed_vol
                    elif len(np.where((dat.loc[Nf_tropha, "p_temp"] < 0) & (dat.loc[Nf_tropha, "qliquid"] < fed_vol))[
                                 0]) > 0:  # if they have less qliquid than fed volume
                        temp = np.unique(np.random.choice(sorted(Nempty), len(
                            (np.where((dat.loc[Nf_tropha, "p_temp"] < 0) & (dat.loc[Nf_tropha, "qliquid"] < fed_vol)))[
                                0])))  # select how many indivs to feed
                        q_feed = dat.loc[Nf_tropha[(dat.loc[Nf_tropha, "p_temp"] < 0) & (
                                dat.loc[Nf_tropha, "qliquid"] < fed_vol)], "qliquid"].values[
                            len(temp) - 1]  # feed them with remaining stuff
                        dat.loc[temp, "fed"] = q_feed / fed_vol  # feed them with remaining stuff
                        dat.loc[temp, "fed"] = np.where((dat.loc[temp, "fed"] > 0.9), 1, dat.loc[
                            temp, "fed"])  # those fed more than 90% are considered full, --- simplified version, even if they have more than 1
                        dat.loc[Nf_tropha, "qliquid"] = np.where(
                            (dat.loc[Nf_tropha, "p_temp"] < 0) & (dat.loc[Nf_tropha, "qliquid"] < fed_vol), 0,
                            dat.loc[Nf_tropha, "qliquid"])  # individuals that passed liquid get to zero
                    dat.loc[:, "p_temp"] = 0  # reboot proba to 0
                    temp = 0
                    dat.loc[Nf_tropha, "qliquid"] = np.where((dat.loc[Nf_tropha, "qliquid"] < 0.1), 0, dat.loc[
                        Nf_tropha, "qliquid"])  # those with less than 10% are considered empty, just in case of miss-adjustment
            ############
            # 2.2 - SB  #
            ############
            if len(Nf_sb) > 0:
                Nempty = np.where((dat.loc[:, "fed"] < 1) & (dat.loc[:, "where"] == 0))[0]
                if (len(Nempty) > 0):  # if empty ants in nest
                    n_feeds = np.random.choice(range(1, 4), len(Nf_sb), replace=True)
                    dat.loc[Nf_sb, "p_temp"] = np.random.random(len(Nf_sb)) - 1 / (2 / (
                            len(Nempty) / N) + 0.524 * n_feeds)  # We multiply prob by 1-3 randomly, but could fix this to X -- can be rewriten as len(Nempty)/(0.524*len(Nempty)*N + 2*n_feeds)
                    q_feed = np.random.random(len(Nf_sb))  # random percentatge of food to pass
                    q_feed = np.where(q_feed > 0.9, 1, q_feed)  # if they pass more than 90%, let's say they pass everything
                    q_feed *= dat.loc[Nf_sb, "qliquid"]  # determine quantity in liquid
                    dat.loc[Nf_sb, "qliquid"] = np.where((dat.loc[Nf_sb, "p_temp"] < 0), dat.loc[Nf_sb, "qliquid"] - q_feed,
                                                         dat.loc[Nf_sb, "qliquid"])  # individuals pass liquid
                    temp = np.sum(n_feeds[np.where(dat.loc[Nf_sb, "p_temp"] < 0)[0]])  # how many individuals are fed
                    if temp <= len(Nempty):
                        temp = np.random.choice(sorted(Nempty), temp,
                                                replace=False)  # choose them randomly among the pool of empty; cannot choose same indiv twice
                        dat.loc[temp, "fed"] += (np.repeat(
                            q_feed.iloc[np.where(dat.loc[Nf_sb, "p_temp"] < 0)[0]].values / n_feeds[
                                np.where(dat.loc[Nf_sb, "p_temp"] < 0)[0]], n_feeds[np.where(dat.loc[Nf_sb, "p_temp"] < 0)[
                                0]])) / fed_vol  # those are fed, we increase fed cell, according to the amount passed and divided by the total numbers of individuals engaged with each forager
                    else:
                        temp = Nempty  # choose them randomly among the pool of empty; cannot choose same indiv twice
                        dat.loc[temp, "fed"] += (np.repeat(
                            q_feed.iloc[np.where(dat.loc[Nf_sb, "p_temp"] < 0)[0]].values / n_feeds[
                                np.where(dat.loc[Nf_sb, "p_temp"] < 0)[0]],
                            n_feeds[np.where(dat.loc[Nf_sb, "p_temp"] < 0)[0]]))[0:len(
                            Nempty)] / fed_vol  # those are fed, we increase fed cell, according to the amount passed and divided by the total numbers of individuals engaged with each forager
                    ###################################
                    #	METHOD SB: COMPLEX VS SIMPLE  #
                    ###################################
                    if method_sb == "complex":
                        dat.loc[temp, "fed"] = np.where((dat.loc[temp, "fed"] > 0.9) & (dat.loc[temp, "fed"] < 1), 1,
                                                        dat.loc[
                                                            temp, "fed"])  # those fed more than 90%, but less than 1 are considered full
                        dat.loc[Nf_sb, "qliquid"] = np.where((dat.loc[Nf_sb, "qliquid"] < 0.1), 0, dat.loc[
                            Nf_sb, "qliquid"])  # those with less than 10% are considered empty
                        if (len(np.where(dat.loc[:, "fed"] > 1)[0]) > 0) & (len(
                                np.where((dat.loc[:, "fed"] < 1) & (dat.loc[:, "where"] == 0))[
                                    0]) > 0):  # if some individuals are overfed and receivers available
                            temp = np.where(dat.loc[:, "fed"] > 1)[0]
                            q_feed = dat.loc[temp, "fed"] - 1
                            dat.loc[temp, "fed"] = 1
                            temp = np.random.choice(
                                sorted(np.where((dat.loc[:, "fed"] < 1) & (dat.loc[:, "where"] == 0))[0]), len(q_feed),
                                replace=True)  # select a random ant to pass food
                            dat.loc[temp, "fed"] += np.sum(q_feed) / len(np.unique(temp))
                            dat.loc[temp, "fed"] = np.where(dat.loc[temp, "fed"] > 1, 1, dat.loc[
                                temp, "fed"])  # if some of them receive more than 1 of food, we set them to 1
                        elif (len(np.where(dat.loc[:, "fed"] > 1)[0]) > 0) & (len(
                                np.where((dat.loc[:, "fed"] < 1) & (dat.loc[:, "where"] == 0))[
                                    0]) == 0):  # if last individuals are overfed
                            dat.loc[dat.loc[:, "fed"] > 1, "fed"] = 1
                    elif method_sb == "simple":
                        dat.loc[temp, "fed"] = np.where((dat.loc[temp, "fed"] > 0.9), 1, dat.loc[
                            temp, "fed"])  # those fed more than 90% are considered full, --- simplified version, even if they have more than 1
                    else:
                        warning_message = "please choose a method_sb: 'complex' or 'simple', default option is 'simple'"
                        return None, warning_message
                    dat.loc[Nf_sb, "qliquid"] = np.where((dat.loc[Nf_sb, "qliquid"] < 0.1), 0, dat.loc[
                        Nf_sb, "qliquid"])  # those with less than 10% are considered empty
                    dat.loc[:, "p_temp"] = 0  # reboot proba to 0
                    temp = 0
            ##################################
            #	Trophallaxis + SB together	 #
            ##################################
            if len(Nf_trsb) > 0:
                # comment: if ant have more than max liquid tropha, they do social bucket
                #			elif liquid is EQUAL or less than max liquid tropha, they do trophallaxis action.
                antsb = np.where(dat.loc[Nf_trsb, "qliquid"] > v_t)
                antt = np.where(dat.loc[Nf_trsb, "qliquid"] <= v_t)
                if len(antsb[0]) > 0:
                    Nempty = np.where((dat.loc[:, "fed"] < 1) & (dat.loc[:, "where"] == 0))[0]
                    if len(Nempty) > 0:  # if empty ants in nest
                        n_feeds = np.random.choice(range(1, 4), len(Nf_trsb[antsb]), replace=True)
                        dat.loc[Nf_trsb[antsb], "p_temp"] = np.random.random(len(Nf_trsb[antsb])) - 1 / (2 / (
                                len(Nempty) / N) + 0.524 * n_feeds)  # We multiply prob by 1-3 randomly, but could fix this to X -- can be rewriten as len(Nempty)/(0.524*len(Nempty)*N + 2*n_feeds)
                        q_feed = np.random.random(len(Nf_trsb[antsb]))  # random percentatge of food to pass
                        q_feed = np.where(q_feed > 0.9, 1,
                                          q_feed)  # if they pass more than 90%, let's say they pass everything
                        q_feed *= dat.loc[Nf_trsb[antsb], "qliquid"]  # determine quantity in liquid
                        dat.loc[Nf_trsb[antsb], "qliquid"] = np.where((dat.loc[Nf_trsb[antsb], "p_temp"] < 0),
                                                                      dat.loc[Nf_trsb[antsb], "qliquid"] - q_feed, dat.loc[
                                                                          Nf_trsb[
                                                                              antsb], "qliquid"])  # individuals pass liquid
                        temp = np.sum(
                            n_feeds[np.where(dat.loc[Nf_trsb[antsb], "p_temp"] < 0)[0]])  # how many individuals are fed
                        if temp <= len(Nempty):
                            temp = np.random.choice(sorted(Nempty), temp,
                                                    replace=False)  # choose them randomly among the pool of empty; cannot choose same indiv twice
                            dat.loc[temp, "fed"] += (np.repeat(
                                q_feed.iloc[np.where(dat.loc[Nf_trsb[antsb], "p_temp"] < 0)[0]].values / n_feeds[
                                    np.where(dat.loc[Nf_trsb[antsb], "p_temp"] < 0)[0]], n_feeds[
                                    np.where(dat.loc[Nf_trsb[antsb], "p_temp"] < 0)[
                                        0]])) / fed_vol  # those are fed, we increase fed cell, according to the amount passed and divided by the total numbers of individuals engaged with each forager
                        else:
                            temp = Nempty  # choose them randomly among the pool of empty; cannot choose same indiv twice
                            dat.loc[temp, "fed"] += (np.repeat(
                                q_feed.iloc[np.where(dat.loc[Nf_trsb[antsb], "p_temp"] < 0)[0]].values / n_feeds[
                                    np.where(dat.loc[Nf_trsb[antsb], "p_temp"] < 0)[0]],
                                n_feeds[np.where(dat.loc[Nf_trsb[antsb], "p_temp"] < 0)[0]]))[0:len(
                                Nempty)] / fed_vol  # those are fed, we increase fed cell, according to the amount passed and divided by the total numbers of individuals engaged with each forager
                        ###################################
                        #	METHOD SB: COMPLEX VS SIMPLE  #
                        ###################################
                        if method_sb == "complex":
                            dat.loc[temp, "fed"] = np.where((dat.loc[temp, "fed"] > 0.9) & (dat.loc[temp, "fed"] < 1), 1,
                                                            dat.loc[
                                                                temp, "fed"])  # those fed more than 90%, but less than 1 are considered full
                            if (len(np.where(dat.loc[:, "fed"] > 1)[0]) > 0) & (len(
                                    np.where((dat.loc[:, "fed"] < 1) & (dat.loc[:, "where"] == 0))[
                                        0]) > 0):  # if some individuals are overfed and receivers available
                                temp = np.where(dat.loc[:, "fed"] > 1)[0]
                                q_feed = dat.loc[temp, "fed"] - 1
                                dat.loc[temp, "fed"] = 1
                                temp = np.random.choice(
                                    sorted(np.where((dat.loc[:, "fed"] < 1) & (dat.loc[:, "where"] == 0))[0]), len(q_feed),
                                    replace=True)  # select a random ant to pass food
                                dat.loc[temp, "fed"] += np.sum(q_feed) / len(np.unique(temp))
                                dat.loc[temp, "fed"] = np.where(dat.loc[temp, "fed"] > 1, 1, dat.loc[
                                    temp, "fed"])  # if some of them receive more than 1 of food, we set them to 1
                            elif (len(np.where(dat.loc[:, "fed"] > 1)[0]) > 0) & (len(
                                    np.where((dat.loc[:, "fed"] < 1) & (dat.loc[:, "where"] == 0))[
                                        0]) == 0):  # if last individuals are overfed
                                dat.loc[dat.loc[:, "fed"] > 1, "fed"] = 1
                        elif method_sb == "simple":
                            dat.loc[temp, "fed"] = np.where((dat.loc[temp, "fed"] > 0.9), 1, dat.loc[
                                temp, "fed"])  # those fed more than 90% are considered full, --- simplified version, even if they have more than 1
                        else:
                            warning_message = "please choose a method_sb: 'complex' or 'simple', default option is 'simple'"
                            return None, warning_message
                        dat.loc[:, "p_temp"] = 0  # reboot proba to 0
                        temp = 0
                        antsb = 0
                    Nempty = np.where((dat.loc[:, "fed"] < 1) & (dat.loc[:, "where"] == 0))[0]
                elif len(antt[0]) > 0:
                    Nempty = np.where((dat.loc[:, "fed"] < 1) & (dat.loc[:, "where"] == 0))[0]
                    if len(Nempty) > 0:  # if empty ants in nest
                        dat.loc[Nf_trsb[antt], "p_temp"] = np.random.random(len(Nf_trsb[antt])) - p_feed_t * len(
                            Nempty) / N  # 1/feeding_time * probability to find empty ant; could be N/(area*something else). Makes up for time passing food --------------- check max_capacity
                        if len(np.where(
                                (dat.loc[Nf_trsb[antt], "p_temp"] < 0) & (dat.loc[Nf_trsb[antt], "qliquid"] >= fed_vol))[
                                   0]) > 0:  # if they have qliquid higher than fed volume
                            temp = random.sample(sorted(Nempty), len((np.where(
                                (dat.loc[Nf_trsb[antt], "p_temp"] < 0) & (dat.loc[Nf_trsb[antt], "qliquid"] >= fed_vol)))[
                                                                         0]))  # we choose randomly which individual is fed with fed_vol
                            dat.loc[temp, "fed"] = 1  # those are fed
                            temp = np.where(
                                (dat.loc[Nf_trsb[antt], "p_temp"] < 0) & (dat.loc[Nf_trsb[antt], "qliquid"] >= fed_vol))
                            dat.loc[Nf_trsb[antt][temp], "p_temp"] = 0
                            dat.loc[Nf_trsb[antt][temp], "qliquid"] = dat.loc[Nf_trsb[antt][temp], "qliquid"] - fed_vol
                            Nempty = np.where((dat.loc[:, "fed"] < 1) & (dat.loc[:, "where"] == 0))[0]
                        elif len(np.where(
                                (dat.loc[Nf_trsb[antt], "p_temp"] < 0) & (dat.loc[Nf_trsb[antt], "qliquid"] < fed_vol))[
                                     0]) > 0:  # if they have less qliquid than fed volume
                            temp = np.unique(np.random.choice(sorted(Nempty), len((np.where(
                                (dat.loc[Nf_trsb[antt], "p_temp"] < 0) & (dat.loc[Nf_trsb[antt], "qliquid"] < fed_vol)))[
                                                                                      0])))  # select how many indivs to feed
                            q_feed = dat.loc[Nf_trsb[antt][(dat.loc[Nf_trsb[antt], "p_temp"] < 0) & (
                                    dat.loc[Nf_trsb[antt], "qliquid"] < fed_vol)], "qliquid"].values[
                                len(temp) - 1]  # feed them with remaining stuff
                            dat.loc[temp, "fed"] = q_feed / fed_vol  # feed them with remaining stuff
                            dat.loc[temp, "fed"] = np.where((dat.loc[temp, "fed"] > 0.9), 1, dat.loc[
                                temp, "fed"])  # those fed more than 90% are considered full, --- simplified version, even if they have more than 1
                            dat.loc[Nf_trsb[antt], "qliquid"] = np.where(
                                (dat.loc[Nf_trsb[antt], "p_temp"] < 0) & (dat.loc[Nf_trsb[antt], "qliquid"] < fed_vol), 0,
                                dat.loc[Nf_trsb[antt], "qliquid"])  # individuals that passed liquid get to zero
                        dat.loc[:, "p_temp"] = 0  # reboot proba to 0
                        temp = 0
                        antt = 0
                dat.loc[Nf_trsb, "qliquid"] = np.where((dat.loc[Nf_trsb, "qliquid"] < 0.1), 0, dat.loc[
                    Nf_trsb, "qliquid"])  # those with less than 10% are considered empty, just in case of miss-adjustment
            # ---------------------------------------------------------#
            # 3- ANTS OUTSIDE phase, exploring arena or going to source#
            # ---------------------------------------------------------#
            # 3.1 Naif ants outside
            if len(Nf_out) > 0:
                # go back to nest?
                dat.loc[Nf_out, "p_temp"] = np.random.random(len(Nf_out)) - p_nest  # probability to go back to nest
                dat.loc[Nf_out, "where"] = np.where((dat.loc[Nf_out, "p_temp"] < 0), 0, 1)  # those that go back to nest
                # find source?
                temp = np.where(dat.loc[Nf_out, "p_temp"] > 0)[0]  # those that are still outside
                Nf_out = Nf_out[temp]
                dat.loc[Nf_out, "p_temp"] = np.random.random(len(Nf_out)) - p_source  # probability to find source
                temp = np.where(dat.loc[Nf_out, "p_temp"] < 0)[0]  # those that get to source
                for i in Nf_out[temp]:
                    if dat.loc[i, "behav"] == 1:
                        dat.loc[i, ["where", "naif", "max_time", "fed"]] = [2, 1, t_t, 1]
                    elif dat.loc[i, "behav"] == 0:
                        dat.loc[i, ["where", "naif", "max_time", "fed"]] = [2, 1, t_sb, 1]
                    else:
                        dat.loc[i, ["where", "naif", "max_time", "fed"]] = [2, 1, t_both, 1]
                dat.loc[:, "p_temp"] = 0
                temp = 0
            # 3.2- for informed
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            # ------- COMMENT: this can be simplified into probability by nº of informed ants foraging, or 1/temps_to_grab or 1/temps distance ----- make another version of the model with proba an compare
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
            temp = np.where(
                (dat.loc[foragers, "timing"] == dat.loc[foragers, "max_time"]) & (dat.loc[foragers, "where"] == 1) & (
                        dat.loc[foragers, "naif"] == 1))[0]  # select those that ARRIVE to source
            for i in temp:
                if dat.loc[i, "behav"] == 1:
                    dat.loc[i, ["where", "timing", "max_time"]] = [2, 0,
                                                                   t_t]  # get into source (2), reboot timing and change max_time for trophallaxis or social bucket, they become fed too
                elif dat.loc[i, "behav"] == 0:
                    dat.loc[i, ["where", "timing", "max_time"]] = [2, 0, t_sb]
                else:
                    dat.loc[i, ["where", "timing", "max_time"]] = [2, 0, t_both]
            dat.loc[foragers, "timing"] = np.where((dat.loc[foragers, "where"] == 1) & (dat.loc[foragers, "naif"] == 1),
                                                   dat.loc[foragers, "timing"] + 1, dat.loc[
                                                       foragers, "timing"])  # who is outside and Informed, increase timing by 1
            dat.loc[:, "p_temp"] = 0  # reboot proba to 0
            temp = 0
            # ---------------------------#
            # 4- ANTS AT THE SOURCE phase#
            # ---------------------------#
            if len(Nf_source) > 0:
                temp = np.where((dat.loc[Nf_source, "timing"] == dat.loc[Nf_source, "max_time"]))[
                    0]  # select those that are full --> do this with probability 1/temps_to_grab or 1/temps distance
                for i in Nf_source[temp]:
                    if dat.loc[i, "behav"] == 1:
                        dat.loc[i, ["where", "timing", "max_time", "qliquid"]] = [3, 0, dist_t,
                                                                                  v_t]  # return (3), reboot timing, max_time & qliquid for trophallaxis or social bucket, use 0 here, no +1 later
                    elif dat.loc[i, "behav"] == 0:
                        dat.loc[i, ["where", "timing", "max_time", "qliquid"]] = [3, 0, dist_sb, v_sb]
                    else:
                        dat.loc[i, ["where", "timing", "max_time", "qliquid"]] = [3, 0, dist_both, v_both]
                dat.loc[np.delete(Nf_source, temp), "timing"] += 1  # who is at source, increase timing by 1
                temp = 0
            # --------------------------------#
            # 5- ANTS GOING BACK TO NEST phase#
            # --------------------------------#
            if len(Nf_ret) > 0:
                temp = np.where((dat.loc[Nf_ret, "timing"] == dat.loc[Nf_ret, "max_time"]))[
                    0]  # select those that get back into the nest --> do this with probability 1/temps_to_grab or 1/temps distance
                for i in Nf_ret[temp]:
                    if dat.loc[i, "behav"] == 1:
                        dat.loc[i, ["where", "timing"]] = [0,
                                                           0]  # return (3), reboot timing, max_time & qliquid for trophallaxis or social bucket, use -1 bc we +1 later
                    elif dat.loc[i, "behav"] == 0:
                        dat.loc[i, ["where", "timing"]] = [0,
                                                           0]  # if we don't change max-time, is one change less when they go out again *** remove this change from top, except if we add randomnes... can be also personality like, each ant always same time, or random choice everytime
                        if random.random() < p_drop_sb * dist_sb:
                            dat.loc[i, "qliquid"] = 0
                        else:
                            pass
                    else:
                        dat.loc[i, ["where", "timing"]] = [0,
                                                           0]  # if we don't change max-time, is one change less when they go out again *** remove this change from top, except if we add randomnes... can be also personality like, each ant always same time, or random choice everytime
                        if random.random() < p_drop_sb * dist_both:
                            dat.loc[i, "qliquid"] -= v_sb
                        else:
                            pass
                dat.loc[np.delete(Nf_ret, temp), "timing"] += 1  # who is returning, increase timing by 1
                temp = 0
            if t in vec_time:
                k += 1
                results[k, 0] = np.sum(dat["fed"])  # number of fed ants
                results[k, 1] = len(np.where(dat.loc[:, "where"] == 0)[0])  # number of ants inside
                results[k, 2] = len(np.where((dat.loc[:, "where"] == 1) | (dat.loc[:, "where"] == 3))[0])  # ants outside
                results[k, 3] = len(np.where(dat.loc[:, "where"] == 2)[0])  # ants at source
                results[k, 4] = len(np.where(dat.loc[:, "naif"] == 1)[0])  # ants that visited the source
                results[k, 5] = t  # timestep
        res = pd.DataFrame(results, columns=["fed", "inside", "outside", "source", "informed", "time", "colony"])
        fin_res.append(res.copy())
    fin_res = pd.concat(fin_res)
    print(time.time() - t0)
    return fin_res, warning_message


''' 
#uncomment this section to run several conditions in a loop

#Choose parameters to test
colonysize=[20, 40, 80, 160]
propforagers=[5, 40]
distancesource=[20, 100]
terraindiff=[0, 1]
sugarcon=[0.1, 0.3, 0.5]
behaviortsb=[0, 1, 2]

import itertools
# Generate all combinations
combi = list(itertools.product(colonysize, propforagers, distancesource, terraindiff, sugarcon, behaviortsb))

columns=['colonysize', 'propforagers', 'distancesource', 'terraindiff', 'sugarcon', 'behaviortsb']
# Create a DataFrame to store the combinations
combi = pd.DataFrame(combi, columns=columns)

# Create empty results DataFrame
columns=['colonysize', 'propforagers', 'distancesource', 'terraindiff', 'sugarcon', 'behaviortsb']
columns += ['fed', 'time', 'colony']
results=pd.DataFrame(columns=columns)

#start iterating
for i in combi.itertuples(index=False):
    res, warn =diacamma_model(N = i.colonysize, Nf = i.propforagers, D= i.distancesource, time_sim = 30, sugar = i.sugarcon, behavior = i.behaviortsb, terrain= i.terraindiff, n_sims = 10)
    #get row with %50 colony fed, or the closest value
    for c in res['colony'].unique():
        temp_res=res.loc[res['colony']==c]
        indt=temp_res.loc[temp_res['fed']<int(i.colonysize/2)].index[-1]
        #if simulation didn't arrive to %50 value we take highest value, to estimate needed time
        if indt == len(temp_res)-1:
            temp=list(temp_res.loc[indt, ['fed', 'time', 'colony']])
        #we take %50 value
        else:
            temp=list(temp_res.loc[indt+1, ['fed', 'time', 'colony']])
        results.loc[len(results)]=[*list(i),*temp] # use * to unpack lists, double [[]] to concat in a single row

results.to_csv('example.csv', index=False)
'''
