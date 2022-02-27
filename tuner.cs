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

public class tuner : MonoBehaviour
{
    AudioClip clip;
    public float samplelength, previoustime = 0f;
    public float[] spectrum = new float[2048];
    float time, truetime;
    List<float> times;
    AudioSource audioSource;
    GameObject error_display;
    void Start()
    {
        times = new List<float>();
        WriteString(null);
        audioSource = GameObject.Find("AudioSource").GetComponent<AudioSource>();
        clip = audioSource.clip;
        error_display = GameObject.Find("error");
        //gameObject.GetComponent<PDA>().main(null, new data(), false);
    }

    void Update()
    {
        clip = audioSource.clip;
        time = audioSource.time;
        truetime += Time.deltaTime;
        float[] samples = new float[getIndexFromTime(2f * samplelength) * clip.channels];

        //AudioListener.GetSpectrumData(spectrum, 0, FFTWindow.Rectangular);
        //spectrum_difference();
        //times.Add(time);

        //for every 0.2 seconds, store the audio data corresponding to the microphone input of the past 0.2 seconds into an array
        if(getIndexFromTime(time-previoustime)>2f*samplelength) clip.GetData(samples, getIndexFromTime(time - 2f * samplelength));
        previoustime = time;
        gameObject.GetComponent<PDA>().main(samples, new data(), false, time);
        Debug.Log(time.ToString()+" note: " + gameObject.GetComponent<PDA>().currentnote);
        if(time != 0f) gameObject.GetComponent<PDA>().overallerror /= time;
        error_display.GetComponent<Text>().text = gameObject.GetComponent<PDA>().overallerror.ToString();
    }

    public int getIndexFromTime(float curTime)
    {
        float lengthPerSample = clip.length / (float)clip.samples;

        return Mathf.FloorToInt(curTime / lengthPerSample);
    }

    void spectrum_difference()
    {
        float average = Queryable.Average(spectrum.AsQueryable());
        List<float> five_largest = new List<float>();
        five_largest = kLargest(spectrum, 5);
        float av_diff = 0f;
        for(int i = 0; i<five_largest.Count(); i++)
        {
            av_diff += five_largest[i] - average;
        }
        av_diff /= five_largest.Count();

        WriteString(time.ToString() + ", "+(av_diff*100).ToString());
    }

    public void WriteString(string text)
    {
        string path = "Assets/test.txt";

        if (text == null) File.Create(path).Close();

        //Write some text to the test.txt file
        StreamWriter writer = new StreamWriter(path, true);

        writer.WriteLine(text);

        writer.Close();
    }

    public static List<float> kLargest(float[] arr, int k)
    {
        Array.Sort(arr);
        Array.Reverse(arr);

        List<float> result = new List<float>();
        for (int i = 0; i < k; i++)
            result.Add(arr[i]);
        return result;
    }
}
