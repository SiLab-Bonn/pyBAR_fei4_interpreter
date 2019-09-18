''' Script to check the correctness of the analysis. The analysis is done on raw data and all results are compared to a recorded analysis.
'''

import os
import unittest
import tables as tb
import numpy as np

from pybar_fei4_interpreter import analysis_utils
from pybar_fei4_interpreter import data_struct
from pybar_fei4_interpreter.data_interpreter import PyDataInterpreter
from pybar_fei4_interpreter.data_histograming import PyDataHistograming


# Get package path
testing_path = os.path.dirname(__file__)  # Get the absoulte path of the online_monitor installation

# Set the converter script path
tests_data_folder = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(testing_path)) + r'/testing/test_analysis_data/'))


def convert_data_array(array, filter_func=None, converter_func=None):  # TODO: add copy parameter, otherwise in-place
    '''Filter and convert raw data numpy array (numpy.ndarray)

    Parameters
    ----------
    array : numpy.array
        Raw data array.
    filter_func : function
        Function that takes array and returns true or false for each item in array.
    converter_func : function
        Function that takes array and returns an array or tuple of arrays.

    Returns
    -------
    data_array : numpy.array
        Data numpy array of specified dimension (converter_func) and content (filter_func)
    '''
    if filter_func:
        array = array[filter_func(array)]
    if converter_func:
        array = converter_func(array)
    return array


def is_data_record(value):
    return np.logical_and(np.logical_and(np.less_equal(np.bitwise_and(value, 0x00FE0000), 0x00A00000), np.less_equal(np.bitwise_and(value, 0x0001FF00), 0x00015000)), np.logical_and(np.not_equal(np.bitwise_and(value, 0x00FE0000), 0x00000000), np.not_equal(np.bitwise_and(value, 0x0001FF00), 0x00000000)))


def get_col_row_array_from_data_record_array(array):
    col, row, _ = get_col_row_tot_array_from_data_record_array(array)
    return col, row


def get_col_row_tot_array_from_data_record_array(array):  # TODO: max ToT
    '''Convert raw data array to column, row, and ToT array

    Parameters
    ----------
    array : numpy.array
        Raw data array.

    Returns
    -------
    Tuple of arrays.
    '''
    def get_col_row_tot_1_array_from_data_record_array(value):
        return np.right_shift(np.bitwise_and(value, 0x00FE0000), 17), np.right_shift(np.bitwise_and(value, 0x0001FF00), 8), np.right_shift(np.bitwise_and(value, 0x000000F0), 4)

    def get_col_row_tot_2_array_from_data_record_array(value):
        return np.right_shift(np.bitwise_and(value, 0x00FE0000), 17), np.add(np.right_shift(np.bitwise_and(value, 0x0001FF00), 8), 1), np.bitwise_and(value, 0x0000000F)

    col_row_tot_1_array = np.column_stack(get_col_row_tot_1_array_from_data_record_array(array))
    col_row_tot_2_array = np.column_stack(get_col_row_tot_2_array_from_data_record_array(array))
    col_row_tot_array = np.vstack((col_row_tot_1_array.T, col_row_tot_2_array.T)).reshape((3, -1), order='F').T  # http://stackoverflow.com/questions/5347065/interweaving-two-numpy-arrays
    try:
        col_row_tot_array_filtered = col_row_tot_array[col_row_tot_array[:, 2] < 14]  # [np.logical_and(col_row_tot_array[:,2]<14, col_row_tot_array[:,1]<=336)]
    except IndexError:
        return np.array([], dtype=np.dtype('>u4')), np.array([], dtype=np.dtype('>u4')), np.array([], dtype=np.dtype('>u4'))
    return col_row_tot_array_filtered[:, 0], col_row_tot_array_filtered[:, 1], col_row_tot_array_filtered[:, 2]  # column, row, ToT


class TestAnalysis(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.interpreter = PyDataInterpreter()
        self.histogram = PyDataHistograming()

    @classmethod
    def tearDownClass(self):  # remove created files
        pass

    def test_libraries_stability(self):  # calls 50 times the constructor and destructor to check the libraries
        for _ in range(50):
            interpreter = PyDataInterpreter()
            histogram = PyDataHistograming()
            del interpreter
            del histogram

    def test_data_alignement(self):  # Test if the data alignment is correct (important to detect data alignment issues)
        hits = np.empty((1,), dtype=[('event_number', np.uint64),
                                     ('trigger_number', np.uint32),
                                     ('trigger_time_stamp', np.uint32),
                                     ('relative_BCID', np.uint8),
                                     ('LVL1ID', np.uint16),
                                     ('column', np.uint8),
                                     ('row', np.uint16),
                                     ('tot', np.uint8),
                                     ('BCID', np.uint16),
                                     ('TDC', np.uint16),
                                     ('TDC_time_stamp', np.uint16),
                                     ('TDC_trigger_distance', np.uint8),
                                     ('trigger_status', np.uint8),
                                     ('service_record', np.uint32),
                                     ('event_status', np.uint16)
                                     ])
        self.assertEqual(self.interpreter.get_hit_size(), hits.itemsize)

    def test_analysis_utils_get_n_cluster_in_events(self):  # check compiled get_n_cluster_in_events function
        event_numbers = np.array([[0, 0, 1, 2, 2, 2, 4, 4000000000, 4000000000, 40000000000, 40000000000], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype=np.int64)  # use data format with non linear memory alignment
        result = analysis_utils.get_n_cluster_in_events(event_numbers[0])
        self.assertListEqual([0, 1, 2, 4, 4000000000, 40000000000], result[:, 0].tolist())
        self.assertListEqual([2, 1, 3, 1, 2, 2], result[:, 1].tolist())

    def test_analysis_utils_get_events_in_both_arrays(self):  # check compiled get_events_in_both_arrays function
        event_numbers = np.array([[0, 0, 2, 2, 2, 4, 5, 5, 6, 7, 7, 7, 8], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype=np.int64)
        event_numbers_2 = np.array([1, 1, 1, 2, 2, 2, 4, 4, 4, 7], dtype=np.int64)
        result = analysis_utils.get_events_in_both_arrays(event_numbers[0], event_numbers_2)
        self.assertListEqual([2, 4, 7], result.tolist())

    def test_analysis_utils_get_max_events_in_both_arrays(self):  # check compiled get_max_events_in_both_arrays function
        event_numbers = np.array([[0, 0, 1, 1, 2], [0, 0, 0, 0, 0]], dtype=np.int64)
        event_numbers_2 = np.array([0, 3, 3, 4], dtype=np.int64)
        result = analysis_utils.get_max_events_in_both_arrays(event_numbers[0], event_numbers_2)
        self.assertListEqual([0, 0, 1, 1, 2, 3, 3, 4], result.tolist())

    def test_map_cluster(self):  # Check the compiled function against result
        cluster = np.zeros((20, ), dtype=tb.dtype_from_descr(data_struct.ClusterInfoTable))
        result = np.zeros((20, ), dtype=tb.dtype_from_descr(data_struct.ClusterInfoTable))
        result[1]["event_number"], result[3]["event_number"], result[4]["event_number"], result[7]["event_number"] = 1, 2, 3, 4

        for index in range(cluster.shape[0]):
            cluster[index]["event_number"] = index

        common_event_number = np.array([0, 1, 1, 2, 3, 3, 3, 4, 4], dtype=np.int64)
        self.assertTrue(np.all(analysis_utils.map_cluster(common_event_number, cluster) == result[:common_event_number.shape[0]]))

    def test_hit_histograming(self):
        raw_data = np.array([67307647, 67645759, 67660079, 67541711, 67718111, 67913663, 67914223, 67847647, 67978655, 68081199, 68219119, 68219487, 68425615, 68311343, 68490719, 68373295, 68553519, 68693039, 68573503, 68709951, 68717058, 68734735, 68604719, 68753999, 68761151, 68847327, 69014799, 69079791, 69211359, 69221055, 69279567, 69499247, 69773183, 69788527, 69998559, 69868559, 69872655, 70003599, 69902527, 70274575, 70321471, 70429983, 70563295, 70574959, 70447631, 70584591, 70783023, 71091999, 70972687, 70985087, 71214815, 71382623, 71609135, 71643519, 71720527, 71897695, 72167199, 72040047, 72264927, 72423983, 77471983, 77602863, 77604383, 77485295, 77616415, 77618927, 77619231, 77639983, 77655871, 77544159, 77548303, 77338399, 77345567, 77346287, 77360399, 77255407, 77386211, 77268287, 77279215, 77409599, 77075983, 76951903, 76980527, 77117023, 76991055, 77011007, 77148127, 77148815, 76827167, 76700031, 76868895, 76758575, 76889567, 76558303, 76429599, 76584783, 76468191, 76610943, 76613743, 76620879, 76629375, 76285999, 76321908, 76194319, 76205599, 76233759, 76065391, 76075839, 76093759, 75801311, 75826319, 75829215, 75699231, 75403887, 75565039, 75439135, 75111711, 75115151, 75251487, 75258399, 75138015, 75303471, 74974111, 74868559, 75030047, 75050079, 74714591, 74722847, 74595103, 74649935, 74656815, 74796511, 74455519, 74391519, 74402607, 74534383, 74189695, 74064911, 74246271, 74116063, 74248719, 74133119, 73935183, 73941087, 73811295, 73663583, 73743423, 73449647, 73453391, 73323743, 73343471, 73474159, 73345087, 73206751, 72899295, 72958559, 72828447, 72542623, 82383232, 67374687, 67503967, 67766575, 68179999, 68052847, 68198239, 68104495, 68235759, 68238223, 68472415, 68490463, 68501279, 68621071, 68623903, 68821791, 68988639, 68864047, 69003183, 68876015, 69007423, 68891407, 69267743, 69272367, 69159567, 69666911, 69684447, 70003247, 70018895, 69898927, 69938543, 69942031, 70198863, 70339919, 70587455, 70462783, 70597679, 70796399, 70800015, 70703887, 71121183, 71323151, 71243535, 71578703, 71467695, 71622879, 71629359, 71831264, 71836511, 71710319, 71992943, 72353855, 72355039, 77606628, 77608287, 77622047, 77510223, 77653263, 77664319, 77546223, 77677471, 77549375, 77213519, 77219551, 77232207, 77234991, 77366511, 77373791, 77389647, 77404383, 77070655, 77087199, 76956975, 76996431, 77009183, 77015327, 76683567, 76840351, 76862255, 76888804, 76548975, 76554767, 76427087, 76560159, 76451967, 76456847, 76468015, 76627295, 76352831, 76354863, 76365887, 75923999, 76074175, 75955439, 76086063, 75774239, 75781535, 75792671, 75662111, 75793647, 75797167, 75827023, 75696543, 75390527, 75522031, 75533663, 75541775, 75432255, 75571535, 75115535, 75247999, 75145197, 75151391, 75160799, 74974991, 74852831, 74871839, 74882783, 75023199, 74896943, 75028767, 75046431, 74922463, 74725711, 74621199, 74658623, 74663183, 74336383, 74484559, 74364526, 74370287, 74370639, 74517983, 74393615, 74205471, 74217359, 74227263, 74231727, 74102559, 74237999, 74248735, 73953599, 73868591, 74000703, 74002975, 73877295, 73664910, 73695967, 73704751, 73579583, 73582639, 73719055, 73405998, 73448207, 73481951, 73008831, 73175087, 73044495, 73058863, 73194895, 73197919, 73093151, 72895567, 72918543, 72947039, 72957919, 82383481, 67392015, 67303135, 67312799, 67318303, 67453727, 67454767, 67634719, 67645887, 67717391, 67914111, 67947919, 67818463, 68052959, 68097215, 68500543, 68711909, 68584735, 68726975, 68741679, 68615471, 68750559, 68755487, 68629311, 68764687, 68765648, 68990175, 69022959, 69023727, 69217327, 69547327, 69665839, 69809983, 69814815, 70006831, 70037807, 70055951, 70068511, 70184031, 70323999, 70334687, 70566095, 70588751, 70723935, 71049695, 70952031, 71084831, 71376863, 71256287, 71611039, 71487727, 71618591, 71623999, 71514239, 71891231, 71897327, 71897663, 72036783, 72391487, 77604975, 77608163, 77621327, 77501983, 77635039, 77646559, 77654671, 77655695, 77546543, 77678383, 77345471, 77224735, 77375519, 77385519, 77393967, 76944399, 76975663, 77114628, 77115231, 77127525, 77142959, 76677423, 76699967, 76722287, 76857647, 76739039, 76883567, 76891615, 76453343, 76584335, 76590623, 76594607, 76600031, 76611167, 76617743, 76622303, 76285999, 76329231, 76335839, 76348175, 76350351, 76356783, 75910383, 75639343, 75787615, 75660079, 75796895, 75797615, 75692559, 75827999, 75833487, 75836479, 75518943, 75568143, 75278943, 75290271, 75297903, 75309391, 75312479, 75315119, 74852223, 74987055, 74858047, 74992943, 74875439, 75008031, 74885407, 75027743, 75055583, 74927839, 74738719, 74629087, 74767391, 74779295, 74789343, 74791247, 74323183, 74454239, 74349455, 74364751, 74516047, 74528559, 74192207, 74201535, 74084367, 74220511, 74109039, 74263263, 74133215, 73807119, 73945313, 73868148, 74001631, 73536815, 73684815, 73711439, 73275407, 73408799, 73052767, 73190975, 73209823, 72788271, 72960607, 72487647, 82383730, 67407151, 67415583, 67322127, 67523871, 67700959, 67583039, 67905375, 67793199, 68159583, 68237791, 68306479, 68492399], np.uint32)
        interpreter = PyDataInterpreter()
        histograming = PyDataHistograming()
        interpreter.set_trig_count(1)
        interpreter.set_warning_output(False)
        histograming.set_no_scan_parameter()
        histograming.create_occupancy_hist(True)
        interpreter.interpret_raw_data(raw_data)
        interpreter.store_event()
        histograming.add_hits(interpreter.get_hits())
        occ_hist_cpp = histograming.get_occupancy()[:, :, 0]
        col_arr, row_arr = convert_data_array(raw_data, filter_func=is_data_record, converter_func=get_col_row_array_from_data_record_array)
        occ_hist_python, _, _ = np.histogram2d(col_arr, row_arr, bins=(80, 336), range=[[1, 80], [1, 336]])
        self.assertTrue(np.all(occ_hist_cpp == occ_hist_python))

    def test_trigger_data_format(self):
        raw_data = np.array([3611295745, 82411778, 82793472, 82411779, 82794496, 82411780, 82795520, 82379013, 82379014, 82379015, 82379016, 67240383, 82379017, 82379018, 82379019, 82379020, 82379021, 82379022, 82379023, 82379024, 82379025,
                             3611361282, 82380701, 82380702, 82380703, 82380704, 82380705, 82380706, 82380707, 67240383, 82380708, 82380709, 82380710, 82380711, 82380712, 82380713, 82380714, 82380715, 82380716,
                             3611426819, 82381368, 82381369, 82381370, 82381371, 82381372, 82381373, 82381374, 67240367, 82381375, 82381376, 82381377, 82381378, 82381379, 82381380, 82381381, 82381382, 82381383,
                             3611492356, 82382035, 82382036, 82382037, 82382038, 82382039, 82382040, 82382041, 67240383, 82382042, 82382043, 82382044, 82382045, 82382046, 82382047, 82382048, 82382049, 82382050,
                             3611557893, 82383726, 82383727, 82383728, 82383729, 82383730, 82383731, 82383732, 67240367, 82383733, 82383734, 82383735, 82383736, 82383737, 82383738, 82383739, 82383740, 82383741],
                            np.uint32)
        raw_data_tlu = np.array([3611295745, 3611361282, 3611426819, 3611492356, 3611557893], np.uint32)
        interpreter = PyDataInterpreter()
        histograming = PyDataHistograming()
        for i in (0, 1, 2):  # 0: trigger data contains trigger number, 1: trigger data contains time stamp, 2: trigger data contains 15bit time stamp and 16bit trigger number
            interpreter.set_trigger_data_format(i)
            interpreter.set_trig_count(16)
            interpreter.set_warning_output(False)
            histograming.set_no_scan_parameter()
            histograming.create_occupancy_hist(True)
            interpreter.interpret_raw_data(raw_data)
            interpreter.store_event()
            histograming.add_hits(interpreter.get_hits())
            hits = interpreter.get_hits()
            if i == 0:
                trigger_number_ref = raw_data_tlu & 0x7FFFFFFF
                trigger_time_stamp_ref = np.zeros_like(raw_data_tlu)
            elif i == 1:
                trigger_number_ref = np.zeros_like(raw_data_tlu)
                trigger_time_stamp_ref = raw_data_tlu & 0x7FFFFFFF
            elif i == 2:
                trigger_number_ref = raw_data_tlu & 0x0000FFFF
                trigger_time_stamp_ref = (raw_data_tlu & 0x7FFF0000) >> 16
            self.assertTrue(np.all(hits["trigger_number"] == trigger_number_ref))
            self.assertTrue(np.all(hits["trigger_time_stamp"] == trigger_time_stamp_ref))

    def test_analysis_utils_in1d_events(self):  # check compiled get_in1d_sorted function
        event_numbers = np.array([[0, 0, 2, 2, 2, 4, 5, 5, 6, 7, 7, 7, 8], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype=np.int64)
        event_numbers_2 = np.array([1, 1, 1, 2, 2, 2, 4, 4, 4, 7], dtype=np.int64)
        result = event_numbers[0][analysis_utils.in1d_events(event_numbers[0], event_numbers_2)]
        self.assertListEqual([2, 2, 2, 4, 7, 7, 7], result.tolist())

    def test_1d_index_histograming(self):  # check compiled hist_2D_index function
        x = np.random.randint(0, 100, 100)
        shape = (100, )
        array_fast = analysis_utils.hist_1d_index(x, shape=shape)
        array = np.histogram(x, bins=shape[0], range=(0, shape[0]))[0]
        shape = (5, )  # shape that is too small for the indices to trigger exception
        exception_ok = False
        try:
            array_fast = analysis_utils.hist_1d_index(x, shape=shape)
        except IndexError:
            exception_ok = True
        except Exception:  # other exception that should not occur
            pass
        self.assertTrue(exception_ok & np.all(array == array_fast))

    def test_2d_index_histograming(self):  # check compiled hist_2D_index function
        x, y = np.random.randint(0, 100, 100), np.random.randint(0, 100, 100)
        shape = (100, 100)
        array_fast = analysis_utils.hist_2d_index(x, y, shape=shape)
        array = np.histogram2d(x, y, bins=shape, range=[[0, shape[0]], [0, shape[1]]])[0]
        shape = (5, 200)  # shape that is too small for the indices to trigger exception
        exception_ok = False
        try:
            array_fast = analysis_utils.hist_2d_index(x, y, shape=shape)
        except IndexError:
            exception_ok = True
        except Exception:  # other exception that should not occur
            pass
        self.assertTrue(exception_ok & np.all(array == array_fast))

    def test_3d_index_histograming(self):  # check compiled hist_3D_index function
        with tb.open_file(os.path.join(tests_data_folder + '/hist_data.h5'), mode="r") as in_file_h5:
            xyz = in_file_h5.root.HistDataXYZ[:]
            x, y, z = xyz[0], xyz[1], xyz[2]
            shape = (100, 100, 100)
            array_fast = analysis_utils.hist_3d_index(x, y, z, shape=shape)
            array = np.histogramdd(np.column_stack((x, y, z)), bins=shape, range=[[0, shape[0] - 1], [0, shape[1] - 1], [0, shape[2] - 1]])[0]
            shape = (50, 200, 200)  # shape that is too small for the indices to trigger exception
            exception_ok = False
            try:
                array_fast = analysis_utils.hist_3d_index(x, y, z, shape=shape)
            except IndexError:
                exception_ok = True
            except Exception:  # other exception that should not occur
                pass
            self.assertTrue(exception_ok & np.all(array == array_fast))


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestAnalysis)
    unittest.TextTestRunner(verbosity=2).run(suite)
