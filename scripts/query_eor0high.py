#!/usr/bin/env python
# conda activate /data/curtin_mwaeor/sw/conda/dev

import pyvo
from astropy.time import Time, TimeDelta
from sys import stderr, argv

# get gpstime of recent observations, 2 weeks
recent = (Time.now() - TimeDelta(14, format="jd")).gps
tap = pyvo.dal.TAPService("http://vo.mwatelescope.org/mwa_asvo/tap")
obs = (
    tap.search(
        f"""
SELECT
    obs_id, starttime_utc, ra_pointing, dec_pointing, channel_numbers_csv,
    good_tiles, dataquality, sun_elevation,
    deleted_flag, gpubox_files_archived,
    freq_res, int_time, obsname
FROM mwa.observation
WHERE CONTAINS(
    POINT('ICRS', ra_pointing, dec_pointing),  -- pointing center
    CIRCLE('ICRS', 60, -27.0, 5)      -- is 5 degrees off eor1
) = 1
AND channel_numbers_csv LIKE '%137%'           -- has channel 137 (175MHz)
AND obs_id > {recent}                          -- recently observed
AND dataquality <= 1                           -- no known issues
AND sun_elevation < 0                          -- sun is not up
AND deleted_flag!='TRUE'                       -- not deleted
AND gpubox_files_archived > 1                  -- data available
-- AND freq_res <= 10                          -- (optional) 10kHz resolution or less
-- AND int_time <= 1                           -- (optional) 1s integration or less
ORDER BY obs_id DESC
"""
    )
    .to_table()
    .to_pandas()
    .dropna(axis=1, how="all")
)

# Add datestr column in YYYYMMDD format
obs['datestr'] = obs['starttime_utc'].str[:10].str.replace('-', '')
obs['grp'] = 'eor0high' + obs['datestr']


path='obs_grp_eor0high_recent.tsv'
obs[['obs_id', 'grp']].to_csv(path, index=False, sep='\t', header=False)
from os.path import realpath
print(realpath(path))

# for each group, create a different file
for grp in obs['grp'].unique():
    path=f'obsids-{grp}.csv'
    obs[obs['grp'] == grp][['obs_id']].to_csv(path, index=False, sep='\t', header=False)
    print(realpath(path))
