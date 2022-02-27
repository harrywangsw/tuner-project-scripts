using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.Linq;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using UnityEditor;
using System.IO;

public class PDA : MonoBehaviour
{
    AudioSource audioSource;
    public AudioClip clipforPDA;
    public double minHz, maxHz;
    public float maxallowedtuningerror, lengthofindex, minnoteduration, overallerror = 0f;
    public int nCandidates, chuncksize, nResolution;
    double[] results, a;
    float[] samples, sampleschunck;
    int nLowPeriodInSamples, nHiPeriodInSamples, i, j;
    float noiseloudness, previousswitchtime = 0f, time = 0f, samplestarttime = 0f;
    public bool in_practice = false;
    public string previousnote, currentnote;
    List<float> pitches = new List<float>();
    public List<bool> tuningerror = new List<bool>();
    public GameObject accent, errorslider;

    //trying to erease everything from test.txt before writing more stuff on there
    void start()
    {
        WriteString(null);
    }

    //takes a chunck of sample from an audio file that has a start time and end time(though I don't really care abou the end time) and does autocorrelation to the samples. The results are written in outputsong and will be saved if saved
    public void main(float[] samplesin, data outputsong, bool save, float truetime)
    {
        //if (transform.GetChild(0).gameObject.GetComponent<check_mistake>()) samplestarttime = transform.GetChild(0).gameObject.GetComponent<check_mistake>().time;
        Debug.Log("sample start time: " + time);
        previousnote = null;
        tuningerror = new List<bool>();
        lengthofindex = getTimeFromIndex(chuncksize / 2);

        //if there isn't a sampe chunck being passed in (in other word, we are not checking mistakes from a microphone recording), we take the sample from the audioclip that should have been loaded and attached to this script
        if (samplesin == null)
        {
            Debug.Log("extracting PDA data from clip");
            samplesin = new float[clipforPDA.samples * clipforPDA.channels];
            clipforPDA.GetData(samplesin, 0);
        }
        samples = convertomono(samplesin);
        //using the loudness of 3rd sample chunck as the noise level
        noiseloudness = getloudness(samples.SubArray(3 * chuncksize, chuncksize), 0f);

        // note that higher frequency means lower period
        nLowPeriodInSamples = hzToPeriodInSamples(maxHz, clipforPDA.frequency);
        nHiPeriodInSamples = hzToPeriodInSamples(minHz, clipforPDA.frequency);
        if (nHiPeriodInSamples <= nLowPeriodInSamples) throw new Exception("Bad range for pitch detection.");
        if (samples.Length < nHiPeriodInSamples) throw new Exception("Not enough samples.");

        //we divide the samples into chuncks of size (chuncksize) and do autocorrelation in PDAfor1chunck, then do autocorrelation starting at the middle of this chunck 
        for (i = 0; i < (samples.Length / (chuncksize/2))-1; i++)
        {
            if (i * chuncksize / 2 + chuncksize >= samples.Length) { Debug.Log("samples total length: " + samples.Length + " i times chuncksize " + (i * chuncksize / 2).ToString()); return; }
            sampleschunck = samples.SubArray(i * chuncksize/2, chuncksize);
            PDAfor1chunck(sampleschunck, outputsong, truetime);
        }

        //for (i = 0; i < outputsong.notes.Count; i++) { Debug.Log(time + "note: " + outputsong.notes[i] + " time: " + outputsong.time[i] + " dynamics count: " + outputsong.dynamics.Count); }
        //removeshortnotes(outputsong);
        //merge notes that are the same...until we get loudness thingy working
        mergesamenotes(outputsong);
        outputsong.notes.Add("last");
        outputsong.time.Add(getTimeFromIndex(samples.Length));

        WriteString("   \n   \n");
        //for (i = 0; i < outputsong.notes.Count; i++)  dynamics(i, outputsong);
        if (save) saveanalysis.SavePlayer(outputsong, clipforPDA.name);
        pitches = new List<float>();
    }
    

    public void PDAfor1chunck(float[] sampleschunck, data outputsong, float truetime)
    {
        time = getTimeFromIndex(i*chuncksize/2);
        float currentloudness = getloudness(sampleschunck, 0f);
        Debug.Log(time + " loudness: " + currentloudness);
        Debug.Log("wdddd: " + identifynote(492f));
        Tuple<float, float>[] res = autocorrelation(sampleschunck, truetime);
        Debug.Log(truetime + "/ " + identifynote(res[0].Item1) + " " + identifynote(res[1].Item1) + " " + identifynote(res[2].Item1) + " " + identifynote(res[3].Item1) + " " + identifynote(res[4].Item1) + " " + identifynote(res[5].Item1) + " " + identifynote(res[6].Item1) + " " + identifynote(res[7].Item1));
        string currentnote_without_octave;
        if (currentloudness >= 0.01f)
        {
            string[] notes = new string[res.Length];
            int i;
            for (i = 0; i< res.Length; i++)
            {
                notes[i] = identifynote((float)res[i].Item1).Split('_')[0];
                //Debug.Log(time + " wtf: " + notes[i]);
            }
            var groups = notes.GroupBy(v => v);
            int maxCount = groups.Max(g => g.Count());
            currentnote_without_octave = groups.First(g => g.Count() == maxCount).Key;

            //choose the ocatave that has the heaviest weight
            float[] check_octave = new float[8];
            for (i = 0; i < res.Length; i++)
            {
                if (notes[i] == currentnote_without_octave) 
                {
                    check_octave[int.Parse(identifynote(res[i].Item1).Split('_')[1])]+=res[i].Item2;
                }
            }
            currentnote = currentnote_without_octave+"_"+ Array.IndexOf(check_octave, check_octave.Max()).ToString();

            float averagefreq = 0f;
            float weight = 0f;
            for(i = 0; i<res.Length; i++)
            {
                if(identifynote(res[i].Item1) == currentnote)
                {
                    averagefreq += res[i].Item1*res[i].Item2;
                    weight+=res[i].Item2;
                }
            }
            averagefreq /= weight;
            if(!float.IsNaN(averagefreq)) checkpitch(averagefreq);
        }
        else
        {
            currentnote = "rest";
            currentloudness = 0f;
        }

        if (currentnote != previousnote)
        {
            //Debug.Log(time + "note: " + currentnote);
            outputsong.notes.Add(currentnote);
            outputsong.time.Add(time);
            previousswitchtime = time;
            GameObject.Find("current_note").GetComponent<Text>().text = currentnote;
        }

        //WriteString(time + ", " + currentloudness);
        outputsong.dynamics.Add(currentloudness);
        previousnote = currentnote;
    }

    //currentindex: get the index from the dynamics list that corresponds to the note switch time 
    //highestindex: the peak of the note measured by how many indexes AFTER currentindex
    public void dynamics(int index, data outputsong)
    {
        if (outputsong.notes[index] == "rest") return;
        int i, j, nextindex, highestindex = 0;
        float riserate, fallrate, average;
        int risedura, falldura;
        int currentindex = Mathf.FloorToInt(outputsong.time[index] / outputsong.time[outputsong.time.Count - 1] * outputsong.dynamics.Count);
        if (index < outputsong.time.Count - 2) nextindex = Mathf.FloorToInt(outputsong.time[index + 1] / outputsong.time[outputsong.time.Count - 1] * outputsong.dynamics.Count);
        else return;
        float highest = 0f, lowestindex, lowest, highesttime;
        float[] indexes = new float[nextindex-currentindex];

        Debug.Log("index: "+index+" current index: " + currentindex + " length: " + nextindex);
        indexes = outputsong.dynamics.ToArray().SubArray(currentindex, nextindex-currentindex);
        highest = indexes.Max();
        highestindex = indexes.ToList().IndexOf(highest);
        highesttime = outputsong.time[index]+lengthofindex * highestindex;
        average = indexes.Average();
        List<double> xValues = new List<double>();

        highest = 0f;
        for(i=1; i<indexes.Length-1; i++)
        {
            if(indexes[i-1]<indexes[i] && indexes[i + 1] < indexes[i] && indexes[i]>highest)
            {
                highest = indexes[i];
                highestindex = i;
                highesttime = outputsong.time[index]+lengthofindex * i;
            }
            xValues.Add(outputsong.time[index]+lengthofindex *i);
        }
        xValues.Add(outputsong.time[index] + lengthofindex * i);

        risedura = highestindex;
        falldura = nextindex - highestindex;
        indexes = outputsong.dynamics.ToArray().SubArray(currentindex, highestindex);
        double rSquared, intercept, slope;
        //LinearRegression(xValues.ToArray().SubArray(0, highestindex), Array.ConvertAll(indexes, x => (double)x), out rSquared, out intercept, out slope);
        WriteString(outputsong.time[index] + ", " + outputsong.dynamics[currentindex]);

        Debug.Log("dynam: " + outputsong.notes[index] + " highesttime: " + highesttime + " note time: " + outputsong.time[index] + " " + average);
    }

    void checkpitch(float pitch)
    {
        int power = Mathf.RoundToInt(Mathf.Log((pitch / 440.0f), 1.059463094359f));
        //Debug.Log(time + " the ratio between perfect and played: " + average / (440f * Mathf.Pow(1.059463094359f, power)));
        float error = 0f;

        error = -1200*Mathf.Log(440f*(float)Math.Pow(1.059463094359f, power)/pitch, 2f);
        errorslider.GetComponent<RectTransform>().anchoredPosition = new Vector2(312f * (error/50f), errorslider.GetComponent<RectTransform>().anchoredPosition.y);
        Debug.Log("WTF "+ error.ToString()+" "+pitch.ToString());
        overallerror += error * getTimeFromIndex(chuncksize / 2);
        GameObject.Find("freq").GetComponent<Text>().text = pitch.ToString();

        //if (Math.Abs(average / (440f * Mathf.Pow(1.059463094359f, power)) - 1f) > maxallowedtuningerror) tuningerror.Add(true);
        //else tuningerror.Add(false);

        //check for vibrato
        //if(pitches.Count>2) autocorrelation(Array.ConvertAll(pitches.ToArray(), x=>(double)x));
    }

    void removeshortnotes(data outputsong)
    {
        int c = 1;
        Debug.Log(time + " compare: " + outputsong.notes.Count + " " + tuningerror.Count);
        while(c > 0) 
        {
            c = 0;
            for (int i = 0; i < outputsong.time.Count - 1; i++)
            {
                if (Math.Abs(outputsong.time[i + 1] - outputsong.time[i]) < minnoteduration)
                {
                    Debug.Log(outputsong.time[i] + " note discredited: " + outputsong.notes[i]);
                    outputsong.notes.RemoveAt(i);
                    outputsong.time.RemoveAt(i);
                    tuningerror.RemoveAt(i);
                    c++;
                }
            }
        }
    }

    void mergesamenotes(data outputsong)
    {
        List<float> newtime = new List<float>();
        List<string> newnotes = new List<string>();
        int i, j;
        newnotes.Add(outputsong.notes[0]);
        newtime.Add(outputsong.time[0]);
        for(i = 1; i<outputsong.notes.Count; i++)
        {
            if (outputsong.notes[i] != outputsong.notes[i - 1]) { newnotes.Add(outputsong.notes[i]); newtime.Add(outputsong.time[i]); }
        }
        outputsong.notes = newnotes;
        outputsong.time = newtime;
    }

    void spectrumcommonmultiple(float[] spectrum, float time)
    {
        int peakcount = 0, c = 0, i, j, k;
        float loudest = 0f, loudestindex = 0f, herztperbin = clipforPDA.frequency/spectrum.Length;
        float[,] peaks = new float[18, 2];
        for (i = 18; i < spectrum.Length; i++)
        {
            if (spectrum[i] > loudest)
            {
                loudest = spectrum[i];
                loudestindex = i * herztperbin;
            }
            for (j = 1; j <= 18; j++)
            {
                if (spectrum[i + j] < spectrum[i] && spectrum[i - j] < spectrum[i])
                {
                    c++;
                }
                else
                {
                    break;
                }
            }
            if (c >= 18)
            {
                peaks[peakcount, 0] = (float)i * herztperbin;
                //peaks[peakcount, 1] = 70f+20f*(float)Math.Log(spectrum[i], 10f);
                peaks[peakcount, 1] = spectrum[i];
                if (peaks[peakcount, 1] < 0f)
                {
                    peaks[peakcount, 1] = 0f;
                }
                //Debug.Log(time + " peak #" + peakcount.ToString() + ": " + peaks[peakcount, 0].ToString() + " loudness " + peaks[peakcount, 1]);

                peakcount++;

                if (peakcount >= 18) break;
            }
            c = 0;
        }


        float[,] refinedpeaks = peaks.OrderBy(p => p[1]);
        float max = 0f, loudness = 0f, harmonics = 0f, foundamental = 0f;
        int maxindex = 0;
        max = 0f;
        //C:\Users\harry\Downloads\

        float[] weightedaverage = new float[Mathf.RoundToInt(loudestindex / 27.5f) + 1];
        float weighing = 0f, weightedsum = 0f, peaksum = 0f, loudnessweighing = 0f, F0canidate;
        int maxharmonics = 1;
        int count = 0;

        for (i = 1; i <= Mathf.RoundToInt(loudestindex / 27.5f); i++)
        {
            F0canidate = loudestindex / i;
            for (j = 2; j < refinedpeaks.GetLength(0); j++)
            {
                harmonics = F0canidate * j;
                for (k = maxindex; k < refinedpeaks.GetLength(0); k++)
                {
                    peaksum += refinedpeaks[k, 0] * refinedpeaks[k, 1];
                    loudnessweighing += refinedpeaks[k, 1];
                    if (Mathf.Abs((refinedpeaks[k, 0] / harmonics) - 1) < 0.08f)
                    {
                        c = j;
                        count++;
                        //Debug.Log(time + " the " + j + "th harmonics of " + identifynote(F0canidate, false) + " is a peak!");
                        if (j <= 4)
                        {
                            weightedsum += refinedpeaks[k, 1];
                        }
                        else
                        {
                            weightedsum += refinedpeaks[k, 1] * weirdweighing(j);
                        }
                        break;
                    }
                }
            }

            for (j = 1; j <= c; j++)
            {
                if (j <= 4)
                {
                    weighing += 1;
                }
                else
                {
                    weighing += weirdweighing(j);
                }
            }
            weightedaverage[i] = weightedsum / weighing;
            if (weightedaverage[i] > max)
            {
                max = weightedaverage[i];
                foundamental = F0canidate;
            }
            if (count > 0)
            {
                //Debug.Log(time + " Foudamental canidate: " + identifynote(F0canidate, false) + "number of harmonics: " + count + " max harmonics: " + c + " weighted average: " + weightedaverage[i]);
            }
            weightedsum = 0f;
            count = 0;
            c = 0;
            weighing = 0f;
            maxharmonics = 1;
        }

        for (j = 1; j < refinedpeaks.GetLength(0); j++)
        {
            harmonics = foundamental * j;
            for (k = 0; k < refinedpeaks.GetLength(0); k++)
            {
                if (Mathf.Abs((refinedpeaks[k, 0] / harmonics) - 1) < 0.08f)
                {
                    c = j;
                    //Debug.Log(time + " while refining, the " + j + "th harmonics of " + identifynote(foundamental, false) + " is a peak!");
                    if (j <= 4)
                    {
                        weightedsum += refinedpeaks[k, 1] * refinedpeaks[k, 0] / j;
                    }
                    else
                    {
                        weightedsum += weirdweighing(j) * refinedpeaks[k, 1] * refinedpeaks[k, 0] / j;
                    }
                    break;
                }
            }
        }
        //Debug.Log(time + " while refining, the Foudamental canidate: " + identifynote(foundamental, false) + " harmonics count: " + c + "weightedsum: " + weightedsum);
        for (j = 1; j <= c; j++)
        {
            harmonics = foundamental * j;
            for (k = 0; k < refinedpeaks.GetLength(0); k++)
            {
                if (Mathf.Abs((refinedpeaks[k, 0] / harmonics) - 1) < 0.08f)
                {
                    if (j <= 4)
                    {
                        weighing += refinedpeaks[k, 1];
                    }
                    else
                    {
                        weighing += refinedpeaks[k, 1] * weirdweighing(j);
                    }
                }
            }
        }
        //foundamental = weightedsum / weighing;

        Debug.Log(time + " base note before refine: " + identifynote((foundamental)) + " " + foundamental);
    }

    Tuple<float, float>[] autocorrelation(float[] samples, float truetime)
    {
        float[] results = new float[nHiPeriodInSamples - nLowPeriodInSamples];
        //WriteString(truetime.ToString() + "\n\n");
        int period;
        //nlowperiod corrsponds to 65 hz, highperiod corresponds to 3000 hz
        for (period = nLowPeriodInSamples; period < nHiPeriodInSamples; period += nResolution)
        {
            float sum = 0;
            float denominatorsum = 0;
            // for each sample, find correlation. (If they are far apart, small)
            for (int i = 0; i < samples.Length - period; i++)
            {
                sum += samples[i] * samples[i + period];
                //denominatorsum += Math.Pow(samples[i], 2) + Math.Pow(samples[i + period], 2);
            }

            float mean = sum / (float)samples.Length;
            results[period - nLowPeriodInSamples] = mean;
            //WriteString(period.ToString() + ", " + (mean*100000).ToString());
        }
        //find peaks
        //List<int> peakindexs = FindPeaks(results.ToList(), 4);

        // find the best indices
        Tuple<int, float>[] bestIndices = findBestCandidates(nCandidates, ref results); //note findBestCandidates modifies parameter
        // convert back to Hz
        Tuple<float, float>[] res = new Tuple<float, float>[nCandidates];
        for (int i = 0; i < nCandidates; i++) res[i] = Tuple.Create(periodInSamplesToHz(bestIndices[i].Item1 + nLowPeriodInSamples, clipforPDA.frequency), bestIndices[i].Item2);
        return res;
    }

    /*double[] newapproch(double[] samples)
    {
        double[] results = new double[nHiPeriodInSamples - nLowPeriodInSamples];
        for (int period = nLowPeriodInSamples; period < nHiPeriodInSamples; period += nResolution)
        {
            double sum = 0;
            double denominatorsum = 0;
            // for each sample, find correlation. (If they are far apart, small)
            for (int i = 0; i < samples.Length - period; i++)
            {
                sum += 2*samples[i] * samples[i + period];
                denominatorsum += Math.Pow(samples[i], 2) + Math.Pow(samples[i + period], 2);
            }

            double mean = sum / denominatorsum;
            results[period - nLowPeriodInSamples] = mean;
        }

        // find the best indices
        int[] bestIndices = findBestCandidates(nCandidates, ref results); //note findBestCandidates modifies parameter
        // convert back to Hz
        double[] res = new double[nCandidates];
        for (int i = 0; i < nCandidates; i++) res[i] = periodInSamplesToHz(bestIndices[i] + nLowPeriodInSamples, clipforPDA.frequency);
        return res;
    }*/

    /*double[] Amdf(double[] samples)
    {
        double[] results = new double[nHiPeriodInSamples - nLowPeriodInSamples];
        for (int period = nLowPeriodInSamples; period < nHiPeriodInSamples; period += nResolution)
        {
            double sum = 0;
            // for each sample, see how close it is to a sample n away. Then sum these.
            for (int i = 0; i < samples.Length - period; i++)
                sum += Math.Abs(samples[i] - samples[i + period]);

            double mean = sum / (double)samples.Length;
            mean *= -1; //somewhat of a hack. We are trying to find the minimum value, but our findBestCandidates finds the max. value.
            results[period - nLowPeriodInSamples] = mean;
        }

        // find the best indices
        int[] bestIndices = findBestCandidates(nCandidates, ref results); //note findBestCandidates modifies parameter
                                                                          // convert back to Hz
        double[] res = new double[nCandidates];
        for (int i = 0; i < nCandidates; i++)
            res[i] = periodInSamplesToHz(bestIndices[i] + nLowPeriodInSamples, clipforPDA.frequency);
        return res;
    }*/

    Tuple<int, float>[] findBestCandidates(int n, ref float[] inputs)
    {
        if (inputs.Length < n) throw new Exception("Length of inputs is not long enough.");
        Tuple<int, float>[] res = new Tuple<int, float>[n]; // will hold indices with the highest amounts.

        for (int c = 0; c < n; c++)
        {
            // find the highest.
            float fBestValue = float.MinValue;
            int nBestIndex = 0;
            for (int i = 0; i < inputs.Length; i++)
                if (inputs[i] > fBestValue) { nBestIndex = i; fBestValue = inputs[i]; }

            // record this highest value
            res[c] = Tuple.Create(nBestIndex, fBestValue);

            // now blank out that index.
            //Debug.Log(time + " note: " + identifynote((float)periodInSamplesToHz(res[c] + nLowPeriodInSamples, clipforPDA.frequency)) + " indice strength: " + fBestValue);
            inputs[nBestIndex] = float.MinValue;
           
        }
        return res;
    }

    static List<int> FindPeaks(List<float> values, int rangeOfPeaks)
    {
        List<int> peaks = new List<int>();

        int checksOnEachSide = rangeOfPeaks / 2;
        for (int i = 0; i < values.Count; i++)
        {
            float current = values[i];
            IEnumerable<float> range = values;
            if (i > checksOnEachSide)
                range = range.Skip(i - checksOnEachSide);
            range = range.Take(rangeOfPeaks);
            if (current == range.Max())
                peaks.Add(i);
        }
        return peaks;
    }

    private static int hzToPeriodInSamples(double hz, int sampleRate)
    {
        return (int)(1 / (hz / (double)sampleRate));
    }
    private static float periodInSamplesToHz(int period, int sampleRate)
    {
        return 1 / (period / (float)sampleRate);
    }

    float[] convertomono(float[] multiChannelSamples)
    {
        float[] preProcessedSamples = new float[multiChannelSamples.Length/clipforPDA.channels];
        int numProcessed = 0;
        float combinedChannelAverage = 0f;
        for (int i = 0; i < multiChannelSamples.Length; i++)
        {
            combinedChannelAverage += multiChannelSamples[i];

            // Each time we have processed all channels samples for a point in time, we will store the average of the channels combined
            if ((i + 1) % clipforPDA.channels == 0)
            {
                preProcessedSamples[numProcessed] = combinedChannelAverage / clipforPDA.channels;
                numProcessed++;
                combinedChannelAverage = 0f;
            }
        }

        Debug.Log("Combine Channels done");
        return preProcessedSamples;
    }

    float weirdweighing(int harmonics)
    {
        if (harmonics <= 4)
        {
            return 1;
        }
        else
        {
            return ((float)Mathf.Log((harmonics * (float)Math.Sqrt((harmonics + 1) / (double)harmonics)), 1.2599210498948731647672106072782f) - (float)Mathf.Log(((harmonics - 1) * (float)Math.Sqrt((double)harmonics / (harmonics - 1))), 1.2599210498948731647672106072782f));
        }
    }

    public float getTimeFromIndex(int index)
    {
        return ((1 / (float)clipforPDA.frequency) * index);
    }

    public string identifynote(float frequency)
    {
        if (frequency == 0f)
        {
            Debug.Log("frequency is zero, bud");
            return null;
        }
        string note;
        double octave;
        int power = Mathf.RoundToInt(Mathf.Log((frequency / 440.0f), 1.059463094359f));


        octave = 4 + Math.Ceiling((double)power / 12);
        //Debug.Log("power: " + power + " frequency: " + frequency + " time: " + time);

        if (power > 0)
        {

            switch (Mathf.Abs(power) % 12)
            {
                case 0:
                    note = "a natural";
                    break;
                case 1:
                    note = "b flat";
                    break;
                case 2:
                    note = "b natural";
                    break;
                case 3:
                    note = "c natural";
                    break;
                case 4:
                    note = "d flat";
                    break;
                case 5:
                    note = "d natural";
                    break;
                case 6:
                    note = "e flat";
                    break;
                case 7:
                    note = "e natural";
                    break;
                case 8:
                    note = "f natural";
                    break;
                case 9:
                    note = "g flat";
                    break;
                case 10:
                    note = "g natural";
                    break;
                case 11:
                    note = "a flat";
                    break;
                default:
                    note = null;
                    return note;
                    break;
            }
        }
        else
        {
            switch (Mathf.Abs(power) % 12)
            {
                case 0:
                    note = "a natural";
                    break;
                case 11:
                    note = "b flat";
                    break;
                case 10:
                    note = "b natural";
                    break;
                case 9:
                    note = "c natural";
                    break;
                case 8:
                    note = "d flat";
                    break;
                case 7:
                    note = "d natural";
                    break;
                case 6:
                    note = "e flat";
                    break;
                case 5:
                    note = "e natural";
                    break;
                case 4:
                    note = "f natural";
                    break;
                case 3:
                    note = "g flat";
                    break;
                case 2:
                    note = "g natural";
                    break;
                case 1:
                    note = "a flat";
                    break;
                default:
                    note = null;
                    return note;
                    break;
            }
        }
        if (note[0] == 'b') octave -= 1;

        return note + "_" + octave.ToString();
    }

    public float getloudness(float[] samples, double freq)
    {
        double square = 0;
        float mean;
        double root = 0;

        // Calculate square
        for (int i = 0; i < samples.Length; i++)
        {
            square += Math.Pow(samples[i], 2);
        }

        // Calculate Mean
        mean = ((float)square / (float)(samples.Length));

        // Calculate Root
        root = Math.Sqrt(mean);

        double f2 = freq * freq;
        double f4 = freq * freq * freq * freq;
        //root*=10 * Math.Log(1.562339d * f4 / ((f2 + 107.65265d * 107.65265d)*(f2 + 737.86223d * 737.86223d))) / Math.Log(10d)+ 10f * Math.Log(2.242881E+16d * f4 / ((f2 + 20.598997d * 20.598997d) * (f2 + 20.598997d * 20.598997d)* (f2 + 12194.22d * 12194.22d) * (f2 + 12194.22d * 12194.22d))) / Math.Log(10d);
        return (float)root;
    }

    public void WriteString(string text)
    {
        string path = "Assets/test.txt";

        if(text == null) File.Create(path).Close();

        //Write some text to the test.txt file
        StreamWriter writer = new StreamWriter(path, true);

        writer.WriteLine(text);

        writer.Close();
    }

    /// <summary>
    /// Fits a line to a collection of (x,y) points.
    /// </summary>
    /// <param name="xVals">The x-axis values.</param>
    /// <param name="yVals">The y-axis values.</param>
    /// <param name="rSquared">The r^2 value of the line.</param>
    /// <param name="yIntercept">The y-intercept value of the line (i.e. y = ax + b, yIntercept is b).</param>
    /// <param name="slope">The slop of the line (i.e. y = ax + b, slope is a).</param>
    public static void LinearRegression(double[] xVals, double[] yVals, out double rSquared, out double yIntercept, out double slope)
    {
        if (xVals.Length != yVals.Length)
        {
            throw new Exception("Input values should be with the same length.");
        }

        double sumOfX = 0;
        double sumOfY = 0;
        double sumOfXSq = 0;
        double sumOfYSq = 0;
        double sumCodeviates = 0;

        for (var i = 0; i < xVals.Length; i++)
        {
            var x = xVals[i];
            var y = yVals[i];
            sumCodeviates += x * y;
            sumOfX += x;
            sumOfY += y;
            sumOfXSq += x * x;
            sumOfYSq += y * y;
        }

        var count = xVals.Length;
        var ssX = sumOfXSq - ((sumOfX * sumOfX) / count);
        var ssY = sumOfYSq - ((sumOfY * sumOfY) / count);

        var rNumerator = (count * sumCodeviates) - (sumOfX * sumOfY);
        var rDenom = (count * sumOfXSq - (sumOfX * sumOfX)) * (count * sumOfYSq - (sumOfY * sumOfY));
        var sCo = sumCodeviates - ((sumOfX * sumOfY) / count);

        var meanX = sumOfX / count;
        var meanY = sumOfY / count;
        var dblR = rNumerator / Math.Sqrt(rDenom);

        rSquared = dblR * dblR;
        yIntercept = meanY - ((sCo / ssX) * meanX);
        slope = sCo / ssX;
    }
}

public static class Extensions
{
    public static T[] SubArray<T>(this T[] array, int offset, int length)
    {
        //Debug.Log("okkk: " + offset + " " + length);
        T[] result = new T[length];
        Array.Copy(array, offset, result, 0, length);
        return result;
    }
}
