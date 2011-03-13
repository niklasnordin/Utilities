//
//  CalculateButton.m
//
//  Created by Niklas Nordin on 2011-01-03.
//  Copyright 2011 nequam. All rights reserved.
//

#import "CalculateButton.h"

@implementation CalculateButton
- (IBAction)now:(id)sender 
{
	NSString *text = [input stringValue];

	NSString *out = @"invalid number";
	//if ([text canBeConvertedToEncoding:@"%f"] == 0)
	{
		double temperature = [text doubleValue];
		double bertil = 2.0*temperature;
		if (bertil > 0)
		{
			out = [NSString stringWithFormat:@"%f",bertil];
		}
		
	}

	[output setString:out];
	NSLog(@"%@",out);
}
@end
