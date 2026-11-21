"""Compare signed arterial changes at fixed and variable residuals (Figure 13)."""
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
LABELS = {"exercise": "Exercise", "dobutamine": "Dobutamine", "nitric_oxide": "Nitric oxide",
          "prostacyclin": "Prostacyclin"}
KEYS = ("p50", "hill", "rq", "hb", "query_p50_shift")
NOMINAL = (26.8, 2.7, .8, 15., 0.)


def read_scores():
    with (HERE/"external_pressure_scores.csv").open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def summarize(rows):
    summary = []
    for condition in LABELS:
        subset = [row for row in rows if row["condition"] == condition]
        keyed = {tuple(float(row[key]) for key in KEYS): row for row in subset}
        if len(keyed) != 324 or len(subset) != 324:
            raise ValueError(f"Expected a complete grid for {condition}")
        nominal = keyed[NOMINAL]
        change = float(nominal["observed_change"])
        if change == 0:
            raise ValueError("The reported mean changes must be nonzero for this classification")
        baseline = float(nominal["observed"]) - change
        direction = math.copysign(1., change)
        agree = opposite = 0
        for row in subset:
            signed = [(float(row[key]) - baseline)*direction for key in ("lower", "upper")]
            agree += min(signed) > 1e-10
            opposite += max(signed) < -1e-10
        summary.append(dict(condition=condition, baseline_mmHg=baseline, observed_change_mmHg=change,
                            fixed_change_mmHg=float(nominal["fixed_change"]),
                            venous_change_mmHg=change+float(nominal["venous_delta_error"]),
                            lower_change_mmHg=float(nominal["lower"])-baseline,
                            upper_change_mmHg=float(nominal["upper"])-baseline,
                            agreement=agree, undetermined=len(subset)-agree-opposite,
                            opposite=opposite, settings=len(subset)))
    for row in summary:
        others = [item for item in summary if item is not row]
        differences = [abs(item["fixed_change_mmHg"]-item["observed_change_mmHg"])
                       - abs(item["venous_change_mmHg"]-item["observed_change_mmHg"]) for item in others]
        row["MAE_difference_omitting_this_condition_mmHg"] = sum(differences)/len(others)
    return summary


def plot(summary):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"svg.fonttype":"none",
                         "svg.hashsalt":"external-direction"})
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.05), gridspec_kw={"width_ratios":[1.3,1]})
    left, right = axes
    colors = ["#18716d", "#d7dbdf", "#ac4c27"]
    left.axvline(0, color="#888888", linewidth=.8)
    for i,row in enumerate(summary):
        y = len(summary)-i-1
        left.plot([row["lower_change_mmHg"],row["upper_change_mmHg"]],[y,y],lw=4,
                  color="#9fb6c4",solid_capstyle="round")
        left.plot(row["observed_change_mmHg"],y+.13,"D",color="#202b38",ms=5)
        left.plot(row["fixed_change_mmHg"],y,"o",color="#126c73",ms=6)
        left.plot(row["venous_change_mmHg"],y-.15,marker="x",color="#a44026",ms=6)
        start = 0
        for label,color in zip(("agreement","undetermined","opposite"),colors):
            count=row[label]
            right.barh(y,count,left=start,height=.52,color=color)
            if count>20:
                right.text(start+count/2,y,str(count),ha="center",va="center",
                           color="white" if label!="undetermined" else "#202b38",fontsize=9)
            elif count:
                right.annotate(str(count),(start+count/2,y),xytext=(5,1),textcoords="offset points",color=color)
            start += count
    left.set_yticks([3,2,1,0],list(LABELS.values()))
    right.set_yticks([])
    for ax in axes:
        ax.set_ylim(-.55,3.55)
        ax.spines[["top","right","left"]].set_visible(False)
        ax.tick_params(axis="y",length=0)
    left.set_title("A   Changes from baseline",loc="left",fontweight="bold",fontsize=10)
    right.set_title("B   Direction across the model class",loc="left",fontweight="bold",fontsize=10)
    left.set_xlabel("Change in arterial oxygen pressure (mmHg)")
    right.set(xlim=(0,337),xticks=[0,108,216,324],xlabel="Number of gas-conversion settings")
    handles=[Line2D([],[],color="#9fb6c4",lw=4,label="Model range"),
             Line2D([],[],linestyle="none",marker="D",color="#202b38",label="Observed"),
             Line2D([],[],linestyle="none",marker="o",color="#126c73",label="Fixed residual"),
             Line2D([],[],linestyle="none",marker="x",color="#a44026",label="Venous change")]
    left.legend(handles=handles,loc="upper center",bbox_to_anchor=(.50,-.30),frameon=False,ncol=2,
                columnspacing=1,fontsize=8)
    right.legend(handles=[Patch(color=c,label=l) for c,l in zip(colors,("Agrees","Unresolved","Opposes"))],
                 loc="upper center",bbox_to_anchor=(.5,-.30),frameon=False,ncol=3,columnspacing=.8,fontsize=8)
    fig.subplots_adjust(left=.13,right=.98,bottom=.33,top=.87,wspace=.23)
    fig.savefig(HERE/"Fig13.svg",metadata={"Date":None})
    fig.savefig(HERE/"Fig13.tif",dpi=600,pil_kwargs={"compression":"tiff_lzw"})
    plt.close(fig)


if __name__ == "__main__":
    summary = summarize(read_scores())
    with (HERE/"external_direction_summary.csv").open("w",encoding="utf-8",newline="") as stream:
        writer = csv.DictWriter(stream,fieldnames=summary[0])
        writer.writeheader()
        writer.writerows(summary)
    plot(summary)

