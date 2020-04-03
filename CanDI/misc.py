import scipy.stats as stats
from scipy.interpolate import interp1d
from CanDI import data

class PercentagesHandler(object):

    """
    Not exactly sure how this one will work yet.
    Definitely will calculate percentage of cell lines where a gene is transcribed/essential as well as
    above a certain quantile.
    """
    pass

class DESeq2(object):
    pass




class BinTranscription:

    """
    This might not end up getting used.
    Could be an okay, but arbitrary way to bin things for DESeq
    """

    def __init__(self, subset, owner, name, norm=False):
        self.subset = self._get_subset(subset, norm)

    def _handle_series(self, values, style, caller):
        knee, elbow = self._bin_transcription(values)

        methods = {"over": lambda x, y: x.ge(y),
                   "under": lambda x, y: x.le(y),
                   "normal": lambda x, y, z: x.gt(y).lt(z)}

        binned = methods[caller](values, knee)

    def _quant_norm(self, array, subset):
        """
        array with samples in the columns and probes across the rows
        """
        arr = array.values
        normed_arr = np.zeros_like(arr)
        ranked = np.argsort(arr, axis=0)
        normed_arr[ranked, np.arange(arr.shape[1])] = np.mean(arr[ranked, np.arange(arr.shape[1])],axis=1)[:,np.newaxis]
        normed_df = pd.DataFrame(normed_arr, index=array.index, columns=array.columns)

        return normed_df

    @staticmethod
    def _bin_transcription(array):

        assert len(array) > 10
        try:
            pdf = stats.gaussian_kde(x).evaluate(x) #probability distribution
        except (ValueError, np.linalg.LinAlgError):
            return
        cdf = np.cumsum(pdf)
        f = interp1d(x, cdf)
        newx = np.linspace(x.min(), x.max(), 17) #smooth curve for KneeLocator
        y = f(newx)
        y = y/y[-1]
#        knee = kneed.KneeLocator(newx, y, curve='concave', direction='increasing').knee
#        elbow = kneed.KneeLocator(newx, 1-y, curve='concave', direction='decreasing').knee
        if self.owner is "gene":
            setattr(self, "knee", knee)
            setattr(self, "elbow", elbow)
        return (knee, elbow)
