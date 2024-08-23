from topone import TopONE
import argparse
import logging

def get_dataframes():
    """
    Extract sequences from FASTA (r, nr, m)
    HDMat
    Fit Transform
    Append HDMat..?
    Append Homologies

    Create dataframe
    """
    # GERMANY
    xbz = 55
    ba_5_2_1 = 641
    ef_1_3 = 10
    for fn in ge_fnames:
        print(fn)
        if "_recom_" in fn:
            ss = targ.get_sample_sequences(fn)
            recom_size = len(ss)
            hd = targ.get_hdmatrices(ss)
            hm = targ.fit_transform(hd, True)
            ge.append(hm)
            ge_hd.append(hd)
        elif "_nonrecom_" in fn:
            ss = targ.get_sample_sequences(fn)
            nonrecom_size = (300 - recom_size) // 2
            # nonrecom_size = 162
            # Get both ends to account for two lineages
            ss = ss[:nonrecom_size] + ss[-nonrecom_size:]
            hd = targ.get_hdmatrices(ss)
            hm = targ.fit_transform(hd, True)
            ge.append(hm)
            ge_hd.append(hd)
        else:
            ss = targ.get_sample_sequences(fn)
            # BA.2
            # ss[0:246]

            # B.1.617.2
            # ss[246:(246+1634)]

            # XBC.1
            # ss[(246+1634):]
            if len(ss[:ba_5_2_1]) <= nonrecom_size:
                nr1 = ss[:ba_5_2_1]
            else:
                nr1 = random.sample(ss[:ba_5_2_1], nonrecom_size)
            
            if len(ss[ba_5_2_1:(ba_5_2_1 + ef_1_3)]) <= nonrecom_size:
                nr2 = ss[ba_5_2_1:(ba_5_2_1 + ef_1_3)]
            else:
                nr2 = random.sample(ss[ba_5_2_1:(ba_5_2_1 + ef_1_3)], nonrecom_size)
            
            r = ss[(ba_5_2_1 + ef_1_3):]
            ss = nr1 + nr2 + r
            
            hd = targ.get_hdmatrices(ss)
            hm = targ.fit_transform(hd, True)
            ge.append(hm)
            ge_hd.append(hd)
    backup_ge = ge
    pass

def main():
    pass


if __name__  == "__main__":
    main()