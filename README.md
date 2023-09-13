Presto-single-pulse-search-explaination-
```
if __name__ == '__main__':
    if (0):
        # The following is for profiling
        import hotshot
        prof = hotshot.Profile("hotshot_edi_stats")
        prof.runcall(main)
        prof.close()
        # To see the results:
        if (0):
            from presto.hotshot import stats
            s = stats.load("hotshot_edi_stats")
            s.sort_stats("time").print_stats()
    else:
        # print(sys.version)
        main()
```
* In the first line we assigning name to the block that if it is executed as a main prgram or module etc .
* Imported `hotshot` - used to profile the code block i.e statistical infos like running time and etc 
* So in this block of the code we are just calling the main function as conditions under `if(0)` won`t be executed

```
def main():
    parser = OptionParser(usage)
    parser.add_option("-x", "--xwin", action="store_true", dest="xwin",
                      default=False, help="Don't make a postscript plot, just use an X-window")
    parser.add_option("-p", "--noplot", action="store_false", dest="makeplot",
                      default=True, help="Look for pulses but do not generate a plot")
    parser.add_option("-m", "--maxwidth", type="float", dest="maxwidth", default=0.0,
                      help="Set the max downsampling in sec (see below for default)")
    parser.add_option("-t", "--threshold", type="float", dest="threshold", default=5.0,
                      help="Set a different threshold SNR (default=5.0)")
    parser.add_option("-s", "--start", type="float", dest="T_start", default=0.0,
                      help="Only plot events occuring after this time (s)")
    parser.add_option("-e", "--end", type="float", dest="T_end", default=1e9,
                      help="Only plot events occuring before this time (s)")
    parser.add_option("-g", "--glob", type="string", dest="globexp", default=None,
                      help="Process the files from this glob expression")
    parser.add_option("-f", "--fast", action="store_true", dest="fast",
                      default=False, help="Use a faster method of de-trending (2x speedup)")
    parser.add_option("-b", "--nobadblocks", action="store_false", dest="badblocks",
                      default=True, help="Don't check for bad-blocks (may save strong pulses)")
    parser.add_option("-d", "--detrendlen", type="int", dest="detrendfact", default=1,
                      help="Chunksize for detrending (pow-of-2 in 1000s)")
    parser.add_option("-z", "--gzip", action="store_true", dest="gzip",
                      default=False, help="gzip the output .singlepulse files")
    (opts, args) = parser.parse_args()
```
* This is the part of main fucntion in which they are defining list of options for users to use by using `optparse` python module  
* End of completion we used `(opts, args) = parser.parse_args()` to tell `optparse` that we`re done
  
`parse_args()` returns two values:
* options - an object containing values for all of your options—e.g. if --file takes a single string argument, then options.file will be the filename supplied by the user, or None if the user did not supply that option
* args, the list of positional arguments leftover after parsing options
* `Action = store`tells to take the next argument and store it to your choosen destination
and the flags true/false  just to put the dest value to true or false`
* `help` - use to give some info about the command 
* `dest` - destination 
```
 if len(args)==0:
        if opts.globexp==None:
            print(full_usage)
            sys.exit(0)
        else:
            args = []
            for globexp in opts.globexp.split():
                args += glob.glob(globexp)
 useffts = True
 dosearch = True
 if opts.xwin:
    pgplot_device = "/XWIN"
 else:
    pgplot_device = ""

```
* In this block of code `opts.globexp` is checking the variable if given any in command line . eg ` --globexp .txt `
* usefft & dosearch are assigned to be true 
* if `opts.xwin` is true then pgplot_device is assignned to "/XWIN" otherwise "" 
```
    fftlen = 8192     # Should be a power-of-two for best speed
    chunklen = 8000   # Must be at least max_downfact less than fftlen
    assert(opts.detrendfact in [1,2,4,8,16,32])
    detrendlen = opts.detrendfact*1000
    if (detrendlen > chunklen):
        chunklen = detrendlen
        fftlen = int(next2_to_n(chunklen))
    blocks_per_chunk = chunklen // detrendlen
    overlap = (fftlen - chunklen) // 2
    worklen = chunklen + 2*overlap  # currently it is fftlen...

    max_downfact = 30
    default_downfacts = [2, 3, 4, 6, 9, 14, 20, 30, 45, 70, 100, 150, 220, 300]

    if args[0].endswith(".singlepulse") or args[0].endswith(".singlepulse.gz"):
        filenmbase = args[0][:args[0].rfind(".singlepulse")]
        dosearch = False
    elif args[0].endswith(".dat"):
        filenmbase = args[0][:args[0].rfind(".dat")]
    else:
        filenmbase = args[0]
```
* `fftlen` = 8192 (for fourier transform)
* `chunklen` for chunks you wana divide the whole data series 
* assert statement is used to check if a condition is True, and if it's not, it raises an error.
* `detrendfact` should be a factor of 2  (checking if detrendfact have the values mentioned in code )
* `dtrendlen` is the lenght over which you decide to remove trend
* dtrendlen , chunklen , fflten all are recalculated considering dtrendfact user input (next2_to_n is a lib imported from presto )
* `blocks per chunk` is no. of blocks in a chunk of length chunklen
* overlap is how much , worklen is defined whereas  max_downfact and default_downfact are constant
* checking arg[0] contains what type of extension and we will save it in `filenmbase` without extension 
```
# Don't do a search, just read results and plot
    if not dosearch:
        info, DMs, candlist, num_v_DMstr = \
              read_singlepulse_files(args, opts.threshold, opts.T_start, opts.T_end)
        orig_N, orig_dt = int(info.N), info.dt
        obstime = orig_N * orig_dt
    else:
        DMs = []
        candlist = []
        num_v_DMstr = {}

```
* From previous code block if we get a file with extension `.singlepulse` then dosearch is set to false as we have seen .
* Deciding to do a search if dosearch false then it will get values of info, DMs, candlist, num_v_DMstr from fucntion `read_singlepulse_files()` 
  by reading single pulse files
* As they are referring to single pulse files funtion first let look into down below code :
```
def read_singlepulse_files(infiles, threshold, T_start, T_end):
    DMs = []
    candlist = []
    num_v_DMstr = {}
    for ii, infile in enumerate(infiles):
        if infile.endswith(".singlepulse") or infile.endswith(".singlepulse.gz"):
            filenmbase = infile[:infile.rfind(".singlepulse")]
        else:
            filenmbase = infile
        info = infodata.infodata(filenmbase+".inf")
        DMstr = "%.2f"%info.DM
        DMs.append(info.DM)
        num_v_DMstr[DMstr] = 0
        if ii==0:
            info0 = info
        if os.stat(infile)[6]:
            try:
                cands = np.loadtxt(infile)
                if len(cands.shape)==1:w
                    cands = np.asarray([cands])
                for cand in cands:
                    if cand[2] < T_start: continue
                    if cand[2] > T_end: break
                    if cand[1] >= threshold:
                        candlist.append(candidate(*cand))
                        num_v_DMstr[DMstr] += 1
            except:  # No candidates in the file
                IndexError
    DMs.sort()
    return info0, DMs, candlist, num_v_DMstr
```

* In this code block we are trying to read single pulse file(main objective ) 
* Load all the info in cands (here doubt diff between cand and cand )
* some condition and if that satisfies will store it as a candidate
*  `num_v_DMstr` is used to count no. of candidates of partiuclar DM 
* `candiate` class is used to store that value in a candidate list ?
```
 # Loop over the input files
        for filenm in args:
            if filenm.endswith(".dat"):
                filenmbase = filenm[:filenm.rfind(".dat")]
            else:
                filenmbase = filenm
            info = infodata.infodata(filenmbase+".inf")
            DMstr = "%.2f"%info.DM
            DMs.append(info.DM)
            N, dt = int(info.N), info.dt
            obstime = N * dt
            # Choose the maximum width to search based on time instead
            # of bins.  This helps prevent increased S/N when the downsampling
            # changes as the DM gets larger.
            if opts.maxwidth > 0.0:
                downfacts = [x for x in default_downfacts if x*dt <= opts.maxwidth]
            else:
                downfacts = [x for x in default_downfacts if x <= max_downfact]
            if len(downfacts) == 0:
                downfacts = [default_downfacts[0]]
            if (filenm == args[0]):
                orig_N = N
                orig_dt = dt
                if useffts:
                    fftd_kerns = make_fftd_kerns(default_downfacts, fftlen)
``` 
* If `filname` ends with `.dat` file then using  `rfind` to find the `.dat` file index and try to save only name of the file in `filenmbase` if 
  not then `filename` is saved as it is in `filenmbase`.
* Filenmbase we got from previous step is added with an extension `.inf` to extract the informatiom using `info.info()` and save it in `info`
* Extract DM info from info.dm and add it in DMstr and then append DM values in `DMs`
* calculating observation time(`obstime`) using info of `dt`(time step size) and `N`(no. of time bing )
* `maxwidth` is either given by you or default
* it checks if `maxwidth`(which you entered) is > than 0 
  > if yes then compare it with each (`default_sample list` values)*`dt`  and if it is < than `maxwidth` than add it to the list of `downnfacts`
  > if no then store all the values which are less than the `max_downfact` which is already defined as 30.
* If `downsample` list is empty then the first value of `default_sample` if filled.
* From first file `N` and `dt` values are stored in variable  `orig_N` & `orig_dt`
* If `useffts` is true then it runs a function `make_fftd_kerns` by sending argument - `default_downfacts,fftlen`
As we encontered this function `make_fft_kerns` so we will discuss that before moving forward to next code part d
```
def make_fftd_kerns(downfacts, fftlen):
    fftd_kerns = []
    for downfact in downfacts:
        kern = np.zeros(fftlen, dtype=np.float32)
        # These offsets produce kernels that give results
        # equal to scipy.signal.convolve
        if downfact % 2:  # Odd number
            kern[:downfact//2+1] += 1.0
            kern[-(downfact//2):] += 1.0
        else:             # Even number
            kern[:downfact//2+1] += 1.0
            if (downfact > 2):
                kern[-(downfact//2-1):] += 1.0
        # The following normalization preserves the
        # RMS=1 characteristic of the data
        fftd_kerns.append(rfft(kern / np.sqrt(downfact), -1))
    return fftd_kerns
```
* `fftd_kerns` for storing fft kernels 
* `For` loop over downfact values and creating an array of 0 of length `fftlen`(value given in main fucntion )
* For each `downfact` making a shape of kernel for each downfact values depending on if it is odd or even also if you look at the shape it is not like normal a box car for eg if you take downfact = 8 then the first half 5 index will be 1 and last 3 index will be 1 so that the total is 8 so that the last ones  will overlap with the next one.
* Then normalising then taking fft of that .
* Now continution of code from the main funtion
  ```
  if info.breaks:
                offregions = list(zip([x[1] for x in info.onoff[:-1]],
                                 [x[0] for x in info.onoff[1:]]))

                # If last break spans to end of file, don't read it in (its just padding)
                if offregions[-1][1] == N - 1:
                    N = offregions[-1][0] + 1
  ```

* This code identifies breaks or gaps in the data, constructs offregions to represent data regions without breaks, and then checks if the last 
  break extends to the end of the file. If it does, the code updates N to exclude the padding at the end of the data. This ensures that only the 
  meaningful data regions are considered in further processing.

* what is info.break ?
* if it is true then it creates a list named `offregion`
* no idea about this ?
  ```
  # Compute the file length in detrendlens
            roundN = N // detrendlen * detrendlen
            numchunks = roundN // chunklen
            # Read in the file
            print('Reading "%s"...'%filenm)
            timeseries = np.fromfile(filenm, dtype=np.float32, count=roundN)
            # Split the timeseries into chunks for detrending
            numblocks = roundN // detrendlen
            timeseries.shape = (numblocks, detrendlen)
            stds = np.zeros(numblocks, dtype=np.float64)
            # de-trend the data one chunk at a time
            print('  De-trending the data and computing statistics...')
            for ii, chunk in enumerate(timeseries):
                if opts.fast:  # use median removal instead of detrending (2x speedup)
                    tmpchunk = chunk.copy()
                    tmpchunk.sort()
                    med = tmpchunk[detrendlen//2]
                    chunk -= med
                    tmpchunk -= med
                else:
                    # The detrend calls are the most expensive in the program
                    timeseries[ii] = scipy.signal.detrend(chunk, type='linear')
                    tmpchunk = timeseries[ii].copy()
                    tmpchunk.sort()
  ```
* The purpose of this starting calculation to find roundN is to handle cases where the data length (N) is not an exact multiple of detrendlen, 
  ensuring that you work with complete chunks of detrendlen and any leftover data at the end is excluded from further analysis.
* `numchunks` is the how many chunks of `chunkslen` can fit in total obs length `(roundN)`
* we reads the binary data from the specified file (filenm) using `np.filefrom` and loads it into the timeseries array
* splitting the file into chunks for detrending by finding `numblock`.
* Rehaping the timeseries into rows = numblocks column = detrendlen to form a 2D array which will contains data as `eg,` Row 1 contains the 
  first 200 data points (from index 0 to 199). Row 2 contains the next 200 data points (from index 200 to 399).
* `stds` - making 1D array of length numblock .
* de-trend the data one chunk at a time .
* As enumerates returns  a tuple after every iteration as index and the actual data so we assigned a for loop ii(index) and chunk(actual data)
* if `opts.fast`(usefaster method of detrending) is true then remove median values from chunks(chunk we`re copying is time data or frequency data )
  > quicker than more complex detrending methods.
*  if `opts.fast` is false then at every chunk we will use `scipy.signal.detrend(data,type'linear')`
 > we doing sorting for standard deviation 
```
# The following gets rid of (hopefully) most of the 
                # outlying values (i.e. power dropouts and single pulses)
                # If you throw out 5% (2.5% at bottom and 2.5% at top)
                # of random gaussian deviates, the measured stdev is ~0.871
                # of the true stdev.  Thus the 1.0/0.871=1.148 correction below.
                # The following is roughly .std() since we already removed the median
                stds[ii] = np.sqrt((tmpchunk[detrendlen//40:-detrendlen//40]**2.0).sum() /
                                    (0.95*detrendlen))
            stds *= 1.148
```
* Here we remove 5% what is tom and bottom ?
* what is this formula o std ?
* why we correcting for stds? 
* 1.148 correction factor addition in stds
* why we diving detrendln by by 40 bcz we want to take only middle portion of data s outliers are removed .

```
# sort the standard deviations and separate those with
            # very low or very high values
            sort_stds = stds.copy()
            sort_stds.sort()
            # identify the differences with the larges values (this
            # will split off the chunks with very low and very high stds
            locut = (sort_stds[1:numblocks//2+1] -
                     sort_stds[:numblocks//2]).argmax() + 1
            hicut = (sort_stds[numblocks//2+1:] -
                     sort_stds[numblocks//2:-1]).argmax() + numblocks//2 - 2
            std_stds = np.std(sort_stds[locut:hicut])
            median_stds = sort_stds[(locut+hicut)//2]
            print("    pseudo-median block standard deviation = %.2f" % (median_stds))
            if (opts.badblocks):
                lo_std = median_stds - 4.0 * std_stds
                hi_std = median_stds + 4.0 * std_stds
                # Determine a list of "bad" chunks.  We will not search these.
                bad_blocks = np.nonzero((stds < lo_std) | (stds > hi_std))[0]
                print("    identified %d bad blocks out of %d (i.e. %.2f%%)" % \
                      (len(bad_blocks), len(stds),
                       100.0*float(len(bad_blocks))/float(len(stds))))
                stds[bad_blocks] = median_stds
            else:
                bad_blocks = []
            print("  Now searching...")

```
* first we copied the stds and sort it in an array `sort_stds` then try to find the max std deviation in the lower and higher part of `sort_stds`
* Then store those indexes of max change and calculate the std of std and store it in std_stds  and also median of this in `median_stds`
* Now searching for bad chunks to not waste time on searching them
* 
```
# Now normalize all of the data and reshape it to 1-D
            timeseries /= stds[:,np.newaxis]
            timeseries.shape = (roundN,)
            # And set the data in the bad blocks to zeros
            # Even though we don't search these parts, it is important
            # because of the overlaps for the convolutions
            for bad_block in bad_blocks:
                loind, hiind = bad_block*detrendlen, (bad_block+1)*detrendlen
                timeseries[loind:hiind] = 0.0
            # Convert to a set for faster lookups below
            bad_blocks = set(bad_blocks)
```
* In this we normalize the time series by diving with their respective `stds` and convert it into 1D
* setting all bad blocks in timeseries to 0
```
# Step through the data
            dm_candlist = []
            for chunknum in range(numnuchunks):
                loind = chunknum*chunklen-overlap
                hiind = (chunknum+1)*chunklen+overlap
                # Take care of beginning and end of file overlap issues
                if (chunknum==0): # Beginning of file
                    chunk = np.zeros(worklen, dtype=np.float32)
                    chunk[overlap:] = timeseries[loind+overlap:hiind]
                elif (chunknum==numchunks-1): # end of the timeseries
                    chunk = np.zeros(worklen, dtype=np.float32)
                    chunk[:-overlap] = timeseries[loind:hiind-overlap]
                else:
                    chunk = timeseries[loind:hiind]
```
* Created  a DM list
* `numchunks` is the no. of chunks in `roundN` of length `chunklen`
* `loin` & `hiind` gives the idex of starting & ending of that chunk  in this overalp concept is not clear ?
* begining and end is dealt differently for copying the data from timeseries to `chunk`
* To understand the expression you can try with the image belowfor better visualisation (i struggled without observation ) 
![overal diagram](https://github.com/Aditya1722/Presto-single-pulse-search-explaination-/assets/73752922/d5652001-e1de-4ff2-aa43-c205da4d31ff)

```
# Make a set with the current block numbers
                lowblock = blocks_per_chunk * chunknum
                currentblocks = set(np.arange(blocks_per_chunk) + lowblock)
                localgoodblocks = np.asarray(list(currentblocks -
                                                   bad_blocks)) - lowblock
                # Search this chunk if it is not all bad
                if len(localgoodblocks):
                    # This is the good part of the data (end effects removed)
                    goodchunk = chunk[overlap:-overlap]
```
* In this code block we dealing with the blocks in a chunk - by finding no. of blocks  i'm away from current chunk from  `lowblock`
* Then for fidning the current chuck index we use the block information as seen above
* `localgoodblock` 
```
# need to pass blocks/chunklen, localgoodblocks
                    # dm_candlist, dt, opts.threshold to cython routine

                    # Search non-downsampled data first
                    # NOTE:  these nonzero() calls are some of the most
                    #        expensive calls in the program.  Best bet would 
                    #        probably be to simply iterate over the goodchunk
                    #        in C and append to the candlist there.
                    hibins = np.flatnonzero(goodchunk>opts.threshold)
                    hivals = goodchunk[hibins]
                    hibins += chunknum * chunklen
                    hiblocks = hibins // detrendlen
                    # Add the candidates (which are sorted by bin)
                    for bin, val, block in zip(hibins, hivals, hiblocks):
                        if block not in bad_blocks:
                            time = bin * dt
                            dm_candlist.append(candidate(info.DM, val, time, bin, 1))
```
* In this block we make a `variable` and store the non zeros values indexes using `np.flatnonzerofrom()` from `goodchunk` with a condition it should be greater than threshold S/N .
* hibins contains the indices of potential candidates within the entire dataset.
* hivals contains the values of these potential candidates.
* hiblocks contains the block numbers where these candidates were found within the current chunk.
* Doesn't bin contains indexes of a chunk ? calculate time of each then store info in `dm_candidate`
*  
```
                    # Prepare our data for the convolution
                    if useffts: fftd_chunk = rfft(chunk, -1)
```
* Doing fft of the chunk where sign tells if you wan  do inverse or fft and stored in `fftd_chunk`
* real to complex
  
```
                    # Now do the downsampling...
                    for ii, downfact in enumerate(downfacts):
                        if useffts: 
                            # Note:  FFT convolution is faster for _all_ downfacts, even 2
                            goodchunk = fft_convolve(fftd_chunk, fftd_kerns[ii],
                                                     overlap, -overlap)
                        else:
                            # The normalization of this kernel keeps the post-smoothing RMS = 1
                            kernel = np.ones(downfact, dtype=np.float32) / \
                                     np.sqrt(downfact)
                            smoothed_chunk = scipy.signal.convolve(chunk, kernel, 1)
                            goodchunk = smoothed_chunk[overlap:-overlap]
                        
                        hibins = np.flatnonzero(goodchunk>opts.threshold)
                        hivals = goodchunk[hibins]
                        hibins += chunknum * chunklen
                        hiblocks = hibins // detrendlen
                        hibins = hibins.tolist()
                        hivals = hivals.tolist()
                        # Now walk through the new candidates and remove those
                        # that are not the highest but are within downfact/2
                        # bins of a higher signal pulse
                        hibins, hivals = prune_related1(hibins, hivals, downfact)
                        # Insert the new candidates into the candlist, but
                        # keep it sorted...
                        for bin, val, block in zip(hibins, hivals, hiblocks):
                            if block not in bad_blocks:
                                time = bin * dt
                                bisect.insort(dm_candlist,
                                              candidate(info.DM, val, time, bin, downfact))
```

* first we are  convolving using the info og  fft  and kernel we have
* cant understand the else part here
* good chunk we get from this will be used to calculate again the hibins, hivals , hiblocks
* prune_related1 function is responsible for removing candidates that are close to other candidates but less significant.
```
# Now walk through the dm_candlist and remove the ones that
            # are within the downsample proximity of a higher
            # signal-to-noise pulse
            dm_candlist = prune_related2(dm_candlist, downfacts)
            print("  Found %d pulse candidates"%len(dm_candlist))
            
            # Get rid of those near padding regions
            if info.breaks: prune_border_cases(dm_candlist, offregions)

            # Write the pulses to an ASCII output file, otherwise make an empty file
            if opts.gzip:
                with gzip.GzipFile(filenmbase+'.singlepulse.gz', mode='w') as outfile:
                    if len(dm_candlist):
                        outfile.write("# DM      Sigma      Time (s)     Sample    Downfact\n".encode('utf-8'))
                        for cand in dm_candlist:
                            outfile.write(str(cand).encode('utf-8'))
            else: # not gzipped
                with open(filenmbase+'.singlepulse', mode='w') as outfile:
                    if len(dm_candlist):
                        outfile.write("# DM      Sigma      Time (s)     Sample    Downfact\n")
                        for cand in dm_candlist:
                            outfile.write(str(cand))

            # Add these candidates to the overall candidate list
            for cand in dm_candlist:
                candlist.append(cand)
            num_v_DMstr[DMstr] = len(dm_candlist)
```
