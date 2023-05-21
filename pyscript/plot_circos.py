import os
import pycircos
import matplotlib as mpl
import matplotlib.pyplot as plt
from Bio import SeqIO, SeqFeature
import pandas as pd
import numpy as np
import copy

def get_gbk(fgbk, sample_id=False):
    records = list(SeqIO.parse(fgbk, "genbank"))
    # only return the first record
    record = records[0]
    if sample_id:
        record.id = sample_id
        record.name = sample_id
    record.features.append(SeqFeature.SeqFeature(SeqFeature.FeatureLocation(0, len(record)), type="source", strand=1))
    return record


def get_cds(record):
    # feature locations store in a list [feature_id, start, end, strand]
    feature_locations = []
    for feature in record.features:
        if feature.type == "CDS" or feature.type == "HIT":
            # get phrog annotation if exists
            if "annotation_PHROGS" in feature.qualifiers:
                phrogs = feature.qualifiers["annotation_PHROGS"]
                phrog_integrase = [phrog for phrog in phrogs if "integrase" in phrog]
                if len(phrog_integrase)>0:
                    phrog_anno = phrog_integrase[0]
                else:
                    phrog_anno = phrogs[0]
                phrog_anno = phrog_anno.split("~")[0].split("_")[1]
            else:
                phrog_anno = "No_PHROG"
            feature_locations.append([feature.qualifiers["locus_tag"][0], phrog_anno, feature.location.start.real, feature.location.end.real])

    # convert to dataframe
    df_cds = pd.DataFrame(feature_locations, columns=["feature_id", "phrog", "start", "end"])
    # width of each feature
    df_cds["width"] = df_cds["end"] - df_cds["start"] + 1
    return df_cds


def read_phrog_annot(fin_phrog):
    df_phrog = pd.read_csv(fin_phrog, sep="\t")
    # convert column "phrog" to str
    df_phrog["phrog"] = df_phrog["phrog"].astype(str)
    return df_phrog


def add_phrog_to_cds(df_cds, df_phrog):
    df_cds = df_cds.merge(df_phrog, on="phrog", how="left")
    # fill NaN in column "color" with #222223
    df_cds["color"] = df_cds["color"].fillna("#c9c9c9")
    # fill NaN in column "category" with "Unknown"
    df_cds["category"] = df_cds["category"].fillna("unknown function")
    return df_cds


def plot_circos(record, df_cds_anno, label):
    tick_interval = 10000
    tick_label_size = 8
    window_size = 1000
    sample_id = record.id
    tick_max = len(record.seq)

    garc = pycircos.Garc(arc_id=sample_id, record=record, interspace=0, linewidth=0, facecolor="black", raxis_range=(0,0), label=label, label_visible=True)
    gcircle = pycircos.Gcircle(figsize=(8,8))
    gcircle.add_garc(garc)
    gcircle.set_garcs(0, 360)

    # add tickplot
    tlabels = [f"{tl / 1000:.0f} K" for tl in range(0, tick_max, tick_interval)]
    gcircle.tickplot(sample_id, raxis_range=(600,607), tickinterval=tick_interval, tickcolor="#7B7D7D", ticklabels=tlabels, ticklabelsize=tick_label_size, ticklabelorientation="horizontal")
    gcircle.tickplot(sample_id, raxis_range=(600,602), tickinterval=int(tick_interval/10), tickcolor="#7B7D7D")

    # Plot source
    feat_source = []
    for feat in garc.record.features:
        if feat.type == "source":
            feat_source.append(feat)
            a = feat
    # gcircle.featureplot(sample_id, source=feat_source, raxis_range=(510,511), facecolor="#36454F")

    #Plot CDS
    feat_CDS_plus  = []
    feat_CDS_minus = []
    for feat in garc.record.features:
        if "CDS" in feat.type:
            if feat.strand == 1:
                feat_CDS_plus.append(feat)
            else:
                feat_CDS_minus.append(feat)
    gcircle.featureplot(sample_id, source=feat_CDS_minus, raxis_range=(480,510), facecolor="#C1A7E2")
    gcircle.featureplot(sample_id, source=feat_CDS_plus,  raxis_range=(450,480), facecolor="#69B3E7")


    # PHROGs
    gcircle.barplot(sample_id, data=[1]*len(df_cds_anno["start"]), positions=df_cds_anno["start"], width=df_cds_anno["width"], raxis_range=(513, 598), facecolor=df_cds_anno["color"])

    return gcircle


def main_plot(fgbk, dirout, bcoat_locus_tag, fin_phrog, label="VC-1"):
    record = get_gbk(fgbk)
    df_cds = get_cds(record)
    df_phrog = read_phrog_annot(fin_phrog)
    tick_max = len(record) + 1000
    df_cds_anno = add_phrog_to_cds(df_cds, df_phrog)

    # if feature_id == "BAF3_00080", set color to #FF0000
    df_cds_anno.loc[df_cds_anno["feature_id"] == bcoat_locus_tag, "color"] = "#FF0000"
    # if feature_id == "BAF3_00080", set category to "BCoAT"
    df_cds_anno.loc[df_cds_anno["feature_id"] == bcoat_locus_tag, "category"] = "BCoAT"

    fig1 = plot_circos(record, df_cds_anno, label)
    # fig1.save("circ_BAM3.pdf", format="pdf")
    fout = os.path.join(dirout, label)
    fig1.save(fout, format="pdf")
    return fig1
