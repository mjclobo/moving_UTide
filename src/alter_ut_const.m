% script to alter ut_constants.mat


const.name = [const.name; 'D1SP'; 'D2SP'];
const.freq = [const.freq; (0.0417807462 * 0.96882); (0.0805114007 * 1.00764)];


save ut_constants.mat const sat shallow ut_version