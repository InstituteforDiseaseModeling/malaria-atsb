import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
mpl.rcParams['pdf.fonttype'] = 42

projectdir = os.path.join('E:/', 'Dropbox (IDM)', 'Malaria Team Folder', 'projects', 'atsb')
plotdir = os.path.join(projectdir, 'sim_plots')


def load_data() :

    # from Jake's calibration to HLC data
    ento_calib_count_sim = [
                            2.1116971448063797,
                            0.9107624124735701,
                            5.10931190103281,
                            21.0493044406172,
                            69.85675311088701,
                            92.443541526798,
                            55.11881077289499,
                            15.674506306646899,
                            4.56673455238333
                        ]
    ento_calib_count_ref = [
                        3.43,
                        2.36,
                        3.71,
                        27.64,
                        67.5,
                        92.21,
                        55.29,
                        19.14,
                        6.14
                    ]
    atsb_calib_count_sim = [
                            3.0469523519279,
                            3.4567237570878397,
                            3.2524569295344903,
                            21.069325804710704,
                            47.315703868868,
                            65.343143582345,
                            30.017979621886205,
                            12.163912057876601,
                            4.4418017887508086
                        ]
    atsb_calib_count_ref = [
                        3.0,
                        2.5,
                        1.64,
                        15.64,
                        53.43,
                        67.93,
                        24.71,
                        7.28,
                        1.57
                    ]

    df = pd.DataFrame( { 'control_ref' : ento_calib_count_ref,
                         'control_sim' : ento_calib_count_sim,
                         'intervention_ref' : atsb_calib_count_ref,
                         'intervention_sim' : atsb_calib_count_sim })
    df['month'] = df.index + 4
    return df


if __name__ == '__main__' :

    df = load_data()

    sns.set_style('whitegrid', {'axes.linewidth' : 0.5})
    fig = plt.figure()
    ax = fig.gca()

    for col in [x for x in df.columns.values if '_' in x] :
        ax.plot(df['month'], df[col], label=col)
    ax.legend()
    ax.set_xlabel('month')
    ax.set_ylabel('mosquito count')

    fig.savefig(os.path.join(plotdir, 'jkurlander_HLC_calibration.png'))
    fig.savefig(os.path.join(plotdir, 'jkurlander_HLC_calibration.pdf'), format='PDF')

    plt.show()