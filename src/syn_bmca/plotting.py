xn = model.xn
en = model.en
x_inds = model.x_inds
e_inds = model.e_inds

plt.rcParams["axes.axisbelow"] = False

# fig, ax_matrix = plt.subplots(ncols=2, nrows=2, figsize=(6.5, 5), sharex="row", sharey="row")
fig, ax_matrix = plt.subplots(ncols=2, nrows=2)
for ax in ax_matrix[1, :].flatten():
    ax.set_aspect("equal")

_ = ax_matrix[0, 0].hist(
    xn.values.flatten(), bins=15, lw=1, edgecolor="w", density=True, facecolor=".4"
)
_ = ax_matrix[0, 1].hist(
    np.log(en.values.flatten()), bins=15, lw=1, edgecolor="w", density=True, facecolor=".4"
)

plot_hpd(ax_matrix[1, 0], xn, trace.posterior["chi_ss"][:, :, :, x_inds])
plot_hpd(ax_matrix[1, 1], np.log(en), trace.posterior["log_en_t"][:, :, :, e_inds])

for ax in ax_matrix[1, :]:
    ax.set_rasterization_zorder(1)

ax_matrix[1, 0].set_xlim([-4, 4])
ax_matrix[1, 1].set_xlim([-4, 4])
ax_matrix[1, 0].set_ylim([-4, 4])
ax_matrix[1, 1].set_ylim([-4, 4])

for ax in ax_matrix[0, :]:
    ax.set_xlim([-4, 4])
    ax.set_xticks([-4, -2, 0, 2, 4])

for ax in ax_matrix[1, :]:
    ax.plot([-4, 4], [-4, 4], "--", color=".3", zorder=4, lw=1.5)
    ax.set_xlim([-4, 4])
    ax.set_ylim([-4, 4])

    ax.set_xticks([-4, -2, 0, 2, 4])
    ax.set_yticks([-4, -2, 0, 2, 4])


# ax_matrix[1, 0].fill_between([-1.5, 1.5], [1.5, 1.5], [-1.5, -1.5], zorder=4, color="k", alpha=0.1)

# ax_matrix[0, 0].set_ylim([0, 1.0])

ax_matrix[0, 0].set_title("Metabolomics", fontsize=13)
ax_matrix[0, 1].set_title("Proteomics", fontsize=13)


ax_matrix[0, 0].text(
    0.5,
    1.0,
    r"$\chi$, n={}".format(xn.shape[0] * xn.shape[1]),
    ha="center",
    va="top",
    transform=ax_matrix[0, 0].transAxes,
)
ax_matrix[0, 1].text(
    0.5,
    1.0,
    r"$\log\; \hat{e}$, n=" + str(en.shape[0] * en.shape[1]),
    ha="center",
    va="top",
    transform=ax_matrix[0, 1].transAxes,
)

ax_matrix[0, 1].set_xlabel("Measured")
ax_matrix[-1, 1].set_xlabel("Measured")
ax_matrix[1, 0].set_ylabel("Predicted")
ax_matrix[0, 0].set_ylabel("Frequency")

sns.despine(offset=2.5, trim=True)


def corrwith(left, right, df=True):
    # demeaned data
    left_tiled = np.repeat(left.values[:, np.newaxis, :], right.shape[0], 1)
    right_tiled = np.repeat(right.values[np.newaxis, :, :], left.shape[0], 0)

    ldem = left_tiled - left_tiled.mean(-1)[:, :, np.newaxis]
    rdem = right_tiled - right_tiled.mean(-1)[:, :, np.newaxis]

    num = (ldem * rdem).sum(-1)

    dom = (left.shape[1] - 1) * left_tiled.std(-1) * right_tiled.std(-1)
    correl = num / dom

    if not df:
        return correl
    else:
        return pd.DataFrame(correl, index=left.index, columns=right.index)
