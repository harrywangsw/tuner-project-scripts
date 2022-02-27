using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class start_mic : MonoBehaviour
{
    public int miclooplength;
    void Start()
    {
        string[] microphones;
        AudioSource audioSource;
        audioSource = gameObject.GetComponent<AudioSource>();
        microphones = Microphone.devices;
        audioSource.clip = Microphone.Start(microphones[0], true, miclooplength, 44100);
        while (!(Microphone.GetPosition(microphones[0]) > 0))
        {
            audioSource.Play();
        }
    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
