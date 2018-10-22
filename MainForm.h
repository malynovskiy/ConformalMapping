#pragma once
#include <fstream>
#include <iomanip>
#include <vector>

namespace ConformalMapping 
{
	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	public ref class MainForm : public System::Windows::Forms::Form
	{
	private:
		double** x;
		double** y;

		double Gamma;
		double Q;
		double D;
		int k;
	
		int xOffset;
		int yOffset;
		int zOffset;
	private:
		System::Windows::Forms::PictureBox^  pictureBox1;

		System::Windows::Forms::Button^  button1;
		System::Windows::Forms::Button^  button2;
		System::Windows::Forms::Label^  label1;
		System::Windows::Forms::Label^  label2;
		System::Windows::Forms::Label^  label3;
		System::Windows::Forms::TextBox^  textBox1;
		System::Windows::Forms::TextBox^  textBox2;
		System::Windows::Forms::TextBox^  textBox3;
		System::Windows::Forms::Label^  label4;
		System::Windows::Forms::TextBox^  textBox4;
		System::Windows::Forms::Button^  button3;

		System::Windows::Forms::Button^  scaleUp;
		System::Windows::Forms::Button^  scaleDown;
		System::Windows::Forms::Button^  moveRight;
		System::Windows::Forms::Button^  moveLeft;
		System::Windows::Forms::Button^  moveUp;
		System::Windows::Forms::Button^  moveDown;

		System::ComponentModel::Container ^components;

	public:
		MainForm();

	protected:
		~MainForm();
		
	private:
		void InitData();
		void InitComponents();

		void InitialMapping();
		void Mapping();

		System::Void drawConformalMapping(Color c);
		System::Void showAddInfo();

		System::Void button1_Click(System::Object^  sender, System::EventArgs^  e);
		System::Void button2_Click(System::Object^  sender, System::EventArgs^  e);
		System::Void button3_Click(System::Object^  sender, System::EventArgs^  e);

		System::Void scaleUp_Click(System::Object^  sender, System::EventArgs^  e);
		System::Void scaleDown_Click(System::Object^  sender, System::EventArgs^  e);
		System::Void moveRight_Click(System::Object^  sender, System::EventArgs^  e);
		System::Void moveLeft_Click(System::Object^  sender, System::EventArgs^  e);
		System::Void moveUp_Click(System::Object^  sender, System::EventArgs^  e);
		System::Void moveDown_Click(System::Object^  sender, System::EventArgs^  e);
	
		void copyBuffers(std::vector<std::vector<double>> &tx, std::vector<std::vector<double>> &ty);
	};
}

